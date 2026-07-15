#!/usr/bin/env python3
"""
validate_predictor.py — Validate the multi-zero alignment predictor.

The predictor says these are top candidates in lsq=1000-2000:
  lsq=1237: γ₁(n=16,Δ=0.018) + γ₁₁(n=60,Δ=0.029) + γ₁₂(n=64,Δ=0.032) + γ₁₆(n=76,Δ=0.019)
  lsq=1388: γ₁₁(n=61,Δ=0.0002) + γ₁₂(n=65,Δ=0.0026) + γ₄(n=35,Δ=0.037)
  lsq=1387: γ₁₁(n=61,Δ=0.006) + γ₁₂(n=65,Δ=0.004) + γ₄(n=35,Δ=0.033)

Scan neighborhood ±15 around predicted minima at N=400 (sparse eigsh).
A deep minimum (40s eigsh) validates the predictor.

Also scan neighborhood of lsq=800 (score=1.099 but labeled normal)
to check if score can produce false positives.
"""
import time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 30000
N      = 400
OUTFILE = '/Users/scones/projects/riemann/results/validate_predictor.out'

def vm_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab

def lam_min(lsq, N):
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lsq)
    idx = np.arange(-N, N + 1, dtype=float)
    ns  = idx

    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]; w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI
    d = (TWO_PI * ns / L)[:, None] - grid[None, :]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    arch = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
    B = (32.0 * L * sh * (L ** 2 - 16.0 * np.pi ** 2 * np.outer(ns, ns)) /
         np.outer(L ** 2 + 16.0 * np.pi ** 2 * ns ** 2,
                  L ** 2 + 16.0 * np.pi ** 2 * ns ** 2))

    if VM:
        lkpp = np.array([float(mp.log(k)) for k, _ in VM])
        wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))

    tau = (arch + B - P); tau = (tau + tau.T) / 2.0
    tau_sp = scipy.sparse.csr_matrix(tau)
    try:
        ev, _ = scipy.sparse.linalg.eigsh(tau_sp, k=1, which='SA',
                                           tol=1e-8, maxiter=5000)
        return float(ev[0])
    except Exception:
        ev = scipy.linalg.eigvalsh(tau)
        return float(ev[0])

def scan_region(lsq_center, width, label):
    print(f"\n--- {label} (center={lsq_center}, ±{width}) ---")
    lsq_range = range(lsq_center - width, lsq_center + width + 1, 2)

    with open(OUTFILE, 'a') as f:
        f.write(f"\n--- {label} (center={lsq_center}) ---\n")

    for lsq in lsq_range:
        t0 = time.time()
        lam0 = lam_min(lsq, N)
        dt = time.time() - t0
        deep = " ← DEEP!" if dt > 20 else ""
        line = f"{lsq:6d} {lam0:14.4e}  ({dt:.0f}s){deep}"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")

def main():
    print(f"Predictor validation scan, N={N}")
    print(f"Predicted deep: lsq=1237 (γ₁+γ₁₁+γ₁₂+γ₁₆), lsq=1387-1388 (γ₄+γ₁₁+γ₁₂)")
    print(f"False positive check: lsq=800 (high score but listed as normal)")

    with open(OUTFILE, 'w') as f:
        f.write(f"validate_predictor N={N} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Predicted: lsq=1237, 1387, 1388. False positive check: lsq=800\n\n")

    # Check lsq=800 (high score, should be normal)
    scan_region(800, 10, "lsq=800 (score=1.099, expected normal)")

    # Check predicted deep minima
    scan_region(1237, 15, "lsq=1237 PREDICTION (γ₁+γ₁₆ aligned, score=0.995)")
    scan_region(1388, 10, "lsq=1388 PREDICTION (γ₁₁+γ₁₂ very tight, score=0.929)")

    print(f"\nDone. Output: {OUTFILE}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
scan_2518_cluster.py — Verify the lsq≈2518-2525 mega-resonance cluster.

At lsq=2518:
- γ₇=40.919: n₇=51.0 (Δ=0.0002!) — near-perfect
- γ₈=43.327: n₈=54.0 (Δ=0.002) — excellent
- Plus 6 more zeros with Δ<0.06

This cluster has 8 simultaneously aligned zeros — the best in range 2000-3000.
If the null-space coverage hypothesis is correct, this should be very deep.

Scan lsq=2505-2535 at step=1 with sparse eigsh.
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
OUTFILE = '/Users/scones/projects/riemann/results/scan_2518_cluster.out'

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

def main():
    print(f"Mega-resonance cluster scan: lsq=2505-2535, N={N}")
    print(f"γ₇(n=51,Δ=0.0002) + γ₈(n=54,Δ=0.002) + 6 more zeros")
    print(f"Predicted: deepest resonance in range 2000-3000")
    print(f"{'lsq':>6} {'lam0':>14} {'dt':>5}  deep?")
    print("-" * 50)

    with open(OUTFILE, 'w') as f:
        f.write(f"scan_2518_cluster N={N} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"γ₇(n=51,Δ=0.0002)+γ₈(n=54,Δ=0.002)\n\n")

    for lsq in range(2505, 2536, 1):
        t0 = time.time()
        lam0 = lam_min(lsq, N)
        dt = time.time() - t0
        deep = " ← DEEP!" if dt > 20 else ""
        line = f"{lsq:6d} {lam0:14.4e}  ({dt:.0f}s){deep}"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")

if __name__ == "__main__":
    main()

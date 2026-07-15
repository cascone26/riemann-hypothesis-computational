#!/usr/bin/env python3
"""
gamma4_scan.py — Scan for γ₄ resonance at lsq=γ₄²=925.7

γ₄ = 30.4249, γ₄² = 925.7

Lattice predictions (n from -3 to +2):
  n=-3: 925.7 × e^{-3π/30.425} = 925.7 × 0.733 = 678.5
  n=-2: 925.7 × e^{-2π/30.425} = 925.7 × 0.815 = 754.4
  n=-1: 925.7 × e^{-π/30.425}  = 925.7 × 0.903 = 836.1
  n=0:  925.7                             = 925.7
  n=+1: 925.7 × e^{π/30.425}   = 925.7 × 1.107 = 1025.4

Scan: lsq=700-980 at N=400 using sparse eigsh.
Watch for 40s convergence time as resonance signature.
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
OUTFILE = '/Users/scones/projects/riemann/results/gamma4_scan.out'

GAMMA = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062]

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

def gamma_lattice_label(lsq):
    """Check if lsq is near a lattice point of any of the first 5 zeros."""
    labels = []
    for k, g in enumerate(GAMMA):
        for n in range(-3, 3):
            pred = g**2 * np.exp(n * np.pi / g)
            if abs(lsq - pred) < 3:
                labels.append(f'γ{k+1} n={n:+d}')
    return ' | '.join(labels) if labels else ''

def main():
    lsq_vals = list(range(700, 985, 5))

    print(f"γ₄ resonance scan: lsq=700-980, N={N}")
    print(f"γ₄=30.425, γ₄²=925.7")
    print(f"Predictions: n=-2 at 754, n=-1 at 836, n=0 at 926")
    print(f"{'lsq':>5} {'lam0':>14} {'dt':>5}  label")
    print("-" * 55)

    with open(OUTFILE, 'w') as f:
        f.write(f"gamma4_scan N={N} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"γ₄=30.425, γ₄²=925.7\n\n")

    for lsq in lsq_vals:
        t0 = time.time()
        lam0 = lam_min(lsq, N)
        dt = time.time() - t0
        label = gamma_lattice_label(lsq)
        deep = " ← DEEP!" if dt > 20 else ""
        line = f"{lsq:5d} {lam0:14.4e}  ({dt:.0f}s){deep}  {label}"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")

if __name__ == "__main__":
    main()

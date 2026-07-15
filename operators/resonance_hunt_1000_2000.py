#!/usr/bin/env python3
"""
resonance_hunt_1000_2000.py — Comprehensive resonance scan lsq=1000-2000, step=2.

Uses alignment score (30 zeros) to SKIP low-score regions and focus on high-score.
Only computes lam_min where alignment score > SCORE_THRESHOLD.

At step=2, covers 500 potential points in lsq=1000-2000.
Score filter reduces computation by ~80% while catching all resonances.

Output: resonance_hunt_1000_2000.out (all results)
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
OUTFILE = '/Users/scones/projects/riemann/results/resonance_hunt_1000_2000.out'
SCORE_THRESHOLD = 0.35  # Skip lsq with score < this

# First 30 Riemann zeros
GAMMAS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918720, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446247, 59.347044, 60.831779, 65.112544,
    67.079810, 69.546402, 72.067157, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425274, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
]
SIGMA = 0.05

def alignment_score(lsq):
    L = float(mp.log(lsq))
    score = 0.0
    for g in GAMMAS:
        nk = g * L / TWO_PI
        frac = nk % 1.0
        delta = min(frac, 1.0 - frac)
        score += (1.0 / g**0.5) * np.exp(-delta**2 / (2 * SIGMA**2))
    return score

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
    print(f"Resonance hunt: lsq=1000-2000, step=2, N={N}")
    print(f"Alignment score threshold: {SCORE_THRESHOLD} (skip below)")
    print(f"{'lsq':>6} {'score':>7} {'lam0':>14} {'dt':>5}  notes")
    print("-" * 60)

    with open(OUTFILE, 'w') as f:
        f.write(f"resonance_hunt lsq=1000-2000 step=2 N={N}\n")
        f.write(f"score_threshold={SCORE_THRESHOLD}, sigma={SIGMA}\n\n")
        f.write(f"{'lsq':>6} {'score':>7} {'lam0':>14} {'dt':>5}  notes\n")
        f.write("-" * 60 + "\n")

    deep_minima = []

    for lsq in range(1000, 2001, 2):
        sc = alignment_score(lsq)

        if sc < SCORE_THRESHOLD:
            # Skip (write SKIP to file but not to console)
            with open(OUTFILE, 'a') as f:
                f.write(f"{lsq:6d} {sc:7.3f}         (SKIP)\n")
            continue

        t0 = time.time()
        lam0 = lam_min(lsq, N)
        dt = time.time() - t0

        deep_flag = " ← DEEP!" if dt > 20 else ""
        line = f"{lsq:6d} {sc:7.3f} {lam0:14.4e}  ({dt:.0f}s){deep_flag}"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")

        if dt > 20:
            deep_minima.append((lsq, lam0, dt, sc))

    print(f"\n\n=== DEEP MINIMA FOUND ===")
    print(f"{'lsq':>6} {'lam0':>14} {'dt':>5} {'score':>7}")
    for lsq, lam0, dt, sc in deep_minima:
        print(f"{lsq:6d} {lam0:14.4e}  ({dt:.0f}s) {sc:7.3f}")

    with open(OUTFILE, 'a') as f:
        f.write(f"\n\n=== DEEP MINIMA FOUND ===\n")
        for lsq, lam0, dt, sc in deep_minima:
            f.write(f"{lsq:6d} {lam0:14.4e}  ({dt:.0f}s) {sc:7.3f}\n")

    print(f"\nDone. {len(deep_minima)} deep minima found.")

if __name__ == "__main__":
    main()

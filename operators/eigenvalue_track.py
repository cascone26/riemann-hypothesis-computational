#!/usr/bin/env python3
"""
eigenvalue_track.py — Track individual tau eigenvalues as lambda_sq increases.

Finds the exact lambda_sq values where eigenvalues cross through 0
(or through threshold 1e-3) as lambda_sq increases from 4 to 22.

These crossing points should reveal if dim=2*(lambda_sq-2) jumps
occur at primes, prime powers, or other special values.

Strategy:
- For N=64 (fast), compute tau at lambda_sq = 4.0 to 22.0 in steps of 0.01
- Track each mode's eigenvalue along the path
- Find where eigenvalues become small (near-null crossings)
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
N      = 64

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/eigenvalue_track.json")


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


def build_tau(lambda_sq, N):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
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
        lkpp = np.array([np.log(k) for k, _ in VM])
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
    tau = arch + B - P
    return (tau + tau.T) / 2.0


def main():
    # Coarse scan: lambda_sq from 4 to 22, step 0.05
    lsq_arr = np.arange(4.0, 22.01, 0.05)
    results = []

    THRESH_SMALL = 5e-2  # eigenvalue threshold to report
    N_TOP = 50           # track this many smallest |ev|

    print(f"Scanning {len(lsq_arr)} points with N={N}, G={G}", flush=True)
    t0 = time.time()

    prev_top = None
    for i, lsq in enumerate(lsq_arr):
        tau = build_tau(lsq, N)
        ev  = scipy.linalg.eigvalsh(tau)
        # Sort by magnitude
        order = np.argsort(np.abs(ev))
        top_ev = ev[order[:N_TOP]]
        n_small_02 = int((np.abs(ev) < 0.02).sum())
        n_small_03 = int((np.abs(ev) < 1e-3).sum())

        row = {
            'lsq': round(float(lsq), 5),
            'n_02': n_small_02,
            'n_03': n_small_03,
            'lam_min': float(ev[0]),
            'top50': [float(x) for x in top_ev],
        }
        results.append(row)

        if i % 20 == 0 or (prev_top is not None and n_small_02 != len(results[-2]['top50'][:N_TOP])):
            dt = time.time() - t0
            print(f"  lsq={lsq:.3f}  n<0.02={n_small_02:3}  n<1e-3={n_small_03:3}  "
                  f"lam_min={ev[0]:.4e}  [{i}/{len(lsq_arr)}, {dt:.0f}s]", flush=True)
        prev_top = top_ev

    # Find crossing points: where n_02 changes
    print("\nCrossing points (n_02 jumps):")
    prev_n = results[0]['n_02']
    crossings = []
    for r in results[1:]:
        if r['n_02'] != prev_n:
            crossing = {
                'lsq': r['lsq'],
                'n_02_before': prev_n,
                'n_02_after': r['n_02'],
                'delta': r['n_02'] - prev_n,
            }
            crossings.append(crossing)
            print(f"  lsq={r['lsq']:.4f}: {prev_n} -> {r['n_02']} (delta={r['n_02']-prev_n})")
            prev_n = r['n_02']

    output = {
        'N': N, 'G': G,
        'scan': results,
        'crossings': crossings,
    }
    with open(OUTFILE, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()

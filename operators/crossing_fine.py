#!/usr/bin/env python3
"""
crossing_fine.py — Fine-resolution eigenvalue crossing scan for k=4..10.

Tests whether the second crossing consistently lands at k + 3/4.
Uses step=0.005 (10× finer than eigenvalue_track.py).

Fast because N=64, G=60000, vectorized eigvalsh.
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
                       "../results/crossing_fine.json")


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


def build_tau_eigvalsh(lambda_sq, N):
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
    tau = (tau + tau.T) / 2.0
    return scipy.linalg.eigvalsh(tau)


def main():
    # Fine scan for each k = 4..10
    # Center on the expected crossing region [k+0.2, k+0.9]
    results = {}

    for k in range(4, 11):
        lsq_range = np.arange(k + 0.15, k + 0.95, 0.005)
        print(f"\nScanning k={k} over {len(lsq_range)} points...", flush=True)

        # Baseline dim at k+0.15 (should be dim = 2*(k-2) modes already)
        ev0 = build_tau_eigvalsh(k + 0.15, N)
        n0 = int((np.abs(ev0) < 0.02).sum())
        crossings = []
        prev_n = n0
        pts = []

        t0 = time.time()
        for lsq in lsq_range:
            ev = build_tau_eigvalsh(lsq, N)
            n = int((np.abs(ev) < 0.02).sum())
            pts.append({'lsq': round(float(lsq), 4), 'n': n, 'lam_min': float(ev[0])})
            if n != prev_n:
                crossings.append({
                    'lsq': round(float(lsq), 4),
                    'n_before': prev_n,
                    'n_after': n,
                    'offset': round(float(lsq) - k, 4),
                })
                print(f"  crossing at lsq={lsq:.4f}: {prev_n}->{n}  offset={lsq-k:.4f}", flush=True)
                prev_n = n
        dt = time.time() - t0

        # For each crossing, bisect to ~0.001 precision
        refined = []
        for cr in crossings:
            lsq_lo = cr['lsq'] - 0.005
            lsq_hi = cr['lsq']
            n_target = cr['n_after']
            # Binary search
            for _ in range(10):
                lsq_mid = (lsq_lo + lsq_hi) / 2
                ev_mid = build_tau_eigvalsh(lsq_mid, N)
                n_mid = int((np.abs(ev_mid) < 0.02).sum())
                if n_mid >= n_target:
                    lsq_hi = lsq_mid
                else:
                    lsq_lo = lsq_mid
            refined.append({
                'lsq_lo': round(lsq_lo, 6),
                'lsq_hi': round(lsq_hi, 6),
                'lsq_mid': round((lsq_lo + lsq_hi) / 2, 6),
                'offset': round((lsq_lo + lsq_hi) / 2 - k, 6),
                'is_3over4': abs((lsq_lo + lsq_hi) / 2 - k - 0.75) < 0.01,
            })
            print(f"  refined: offset={refined[-1]['offset']:.6f}  is_3/4={refined[-1]['is_3over4']}", flush=True)

        results[k] = {
            'crossings_coarse': crossings,
            'crossings_refined': refined,
            'scan_dt': round(dt, 1),
        }
        print(f"  k={k}: {len(crossings)} crossings in {dt:.1f}s")

    print("\n\nSummary of second crossings (expected ~k+3/4):")
    for k, v in results.items():
        rc = v['crossings_refined']
        if len(rc) >= 2:
            c2 = rc[1]
            print(f"  k={k}: second crossing at lsq={c2['lsq_mid']:.6f}  offset={c2['offset']:.6f}  "
                  f"offset - 3/4 = {c2['offset'] - 0.75:+.6f}  is_3/4={c2['is_3over4']}")

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'G': G, 'results': {str(k): v for k, v in results.items()}}, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()

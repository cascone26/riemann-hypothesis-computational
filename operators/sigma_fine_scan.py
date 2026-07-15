#!/usr/bin/env python3
"""
sigma_fine_scan.py — Fine scan of sub_cancel vs sigma near 0.5.

Question: Is the sub_cancel minimum at σ=1/2 a SHARP singularity or
a broad minimum? How fast does sub_cancel grow as we move away from σ=0.5?

Tests sigma = 0.40, 0.44, 0.46, 0.48, 0.49, 0.495, 0.499, 0.500, 0.501,
             0.505, 0.51, 0.52, 0.54, 0.56, 0.60 at a fixed lsq=50 (N=150).

At lsq=50, we have enough modes to see the power-law behavior clearly.
"""
import os, json, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
THRESH = 0.02

OUTFILE = os.path.expanduser("~/riemann/results/sigma_fine_scan.json")


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


def build_and_measure(lambda_sq, N_use, sigma):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)

    idx = np.arange(-N_use, N_use + 1, dtype=float)
    ns  = idx
    R   = TWO_PI * N_use / L + 80.0
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
    AB = arch + B
    if VM:
        lkpp = np.array([np.log(k) for k, _ in VM])
        wts  = np.array([lk * (k ** (-sigma)) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))
    tau = AB - P
    tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    null_mask = np.abs(ev) < THRESH
    dim = null_mask.sum()
    if dim == 0:
        return {'dim': 0, 'lam_min': float(ev[0]), 'sub_cancel': float('nan'), 'max_null': float('nan')}
    V_null = V[:, null_mask]
    AB_proj = V_null.T @ AB @ V_null
    P_proj  = V_null.T @ P  @ V_null
    sc = np.linalg.norm(AB_proj - P_proj, 'fro') / (np.linalg.norm(AB_proj, 'fro') + 1e-100)
    return {
        'dim': int(dim),
        'lam_min': float(ev[0]),
        'sub_cancel': float(sc),
        'max_null': float(ev[null_mask].max()),
    }


def main():
    lsq = 50
    N_use = 150  # 3 * lsq
    L = float(mp.log(lsq))

    # Fine scan of sigma
    sigma_vals = [
        0.35, 0.40, 0.42, 0.44, 0.46, 0.48,
        0.490, 0.495, 0.498, 0.499,
        0.500,
        0.501, 0.502, 0.505, 0.510,
        0.52, 0.54, 0.56, 0.58, 0.60, 0.65,
    ]

    print(f"lsq={lsq}, N={N_use}, L={L:.4f}", flush=True)
    print(f"{'sigma':>7}  {'dim':>4}  {'sub_cancel':>12}  {'lam_min':>12}  {'max_null':>10}")

    results = []
    t0 = time.time()
    for sigma in sigma_vals:
        t1 = time.time()
        r = build_and_measure(lsq, N_use, sigma)
        dt = time.time() - t1
        sc = r['sub_cancel']
        print(f"{sigma:7.4f}  {r['dim']:>4}  {sc:>12.4e}  {r['lam_min']:>12.4e}  {r.get('max_null', float('nan')):>10.4e}  ({dt:.1f}s)", flush=True)
        results.append({'sigma': sigma, **r})

    # Find minimum
    sc_vals = [r['sub_cancel'] for r in results if not np.isnan(r['sub_cancel'])]
    sigma_at_min = results[[r['sub_cancel'] for r in results].index(min(sc_vals))]['sigma']
    print(f"\nMinimum sub_cancel at sigma={sigma_at_min:.4f} (sc={min(sc_vals):.4e})")

    # How sharp is the minimum? Measure d(sub_cancel)/d(sigma) on each side of 0.5
    # Compare sc at sigma=0.499 vs 0.500 vs 0.501
    sc_499 = next((r['sub_cancel'] for r in results if r['sigma'] == 0.499), None)
    sc_500 = next((r['sub_cancel'] for r in results if r['sigma'] == 0.500), None)
    sc_501 = next((r['sub_cancel'] for r in results if r['sigma'] == 0.501), None)
    if sc_499 and sc_500 and sc_501:
        print(f"sc(0.499)={sc_499:.4e}, sc(0.500)={sc_500:.4e}, sc(0.501)={sc_501:.4e}")
        ratio_left  = sc_499 / sc_500 if sc_500 > 0 else float('inf')
        ratio_right = sc_501 / sc_500 if sc_500 > 0 else float('inf')
        print(f"sc(0.499)/sc(0.500) = {ratio_left:.2f}")
        print(f"sc(0.501)/sc(0.500) = {ratio_right:.2f}")

    out = {
        'lsq': lsq, 'N': N_use, 'L': round(L, 4), 'thresh': THRESH,
        'results': results,
        'sigma_at_min': sigma_at_min,
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

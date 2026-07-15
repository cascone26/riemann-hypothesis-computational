#!/usr/bin/env python3
"""
dim_scaling_large.py — Verify dim formula at large lsq with proper N = ceil(3*lsq).

Key finding from dim_scaling.py: need N >= lsq (approximately). Here we use N=3*lsq
to ensure the formula holds, then study max_null and sub_cancel power laws.

Runs in ~2-5 min on Olympus.
"""
import os, json, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
from scipy.optimize import curve_fit
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
THRESH = 0.02

OUTFILE = os.path.expanduser("~/riemann/results/dim_scaling_large.json")


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


def build_tau_eigh(lambda_sq, N_use):
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
    ev, V = scipy.linalg.eigh(tau)
    return ev, V, idx, L, arch + B, P


def main():
    # Use N = 3 * lsq to ensure full formula holds
    # lsq values from small to large
    lsq_list = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100]

    results = []
    t0 = time.time()

    for lsq in lsq_list:
        N_use = int(np.ceil(3.0 * lsq))
        expected_dim_plus1 = 2 * (int(lsq) - 2) + 1
        expected_dim_no1   = 2 * (int(lsq) - 2)

        print(f"\nlsq={lsq}, N={N_use} (size {2*N_use+1}×{2*N_use+1})...", flush=True)
        t1 = time.time()
        ev, V, idx, L, AB, P = build_tau_eigh(lsq, N_use)
        dt = time.time() - t1

        null_mask = np.abs(ev) < THRESH
        dim = null_mask.sum()
        lam_min = ev[0]
        max_null = float(ev[null_mask].max()) if dim > 0 else float('nan')

        # Sub_cancel
        V_null = V[:, null_mask]
        AB_proj = V_null.T @ AB @ V_null
        P_proj  = V_null.T @ P  @ V_null
        sc = np.linalg.norm(AB_proj - P_proj, 'fro') / (np.linalg.norm(AB_proj, 'fro') + 1e-100)

        formula_match = "OK+1" if dim == expected_dim_plus1 else ("OK" if dim == expected_dim_no1 else f"GOT {dim}")
        print(f"  dim={dim} expected~{expected_dim_plus1}  [{formula_match}]")
        print(f"  lam_min={lam_min:.4e}  max_null={max_null:.4e}")
        print(f"  sub_cancel={sc:.4e}  L={L:.4f}  ({dt:.1f}s)")

        results.append({
            'lsq': lsq, 'N': N_use, 'L': round(L, 6),
            'dim': int(dim),
            'expected_plus1': expected_dim_plus1,
            'expected_no1': expected_dim_no1,
            'formula_ok': formula_match,
            'lam_min': float(lam_min),
            'max_null': float(max_null),
            'sub_cancel': float(sc),
        })

    # Power law fits
    Ls = np.array([r['L'] for r in results])
    maxnulls = np.array([r['max_null'] for r in results])
    subcancels = np.array([r['sub_cancel'] for r in results])

    fit_mn = fit_sc = None
    try:
        popt, _ = curve_fit(lambda x, C, a: C / x**a, Ls, maxnulls, p0=[0.05, 1.2])
        fit_mn = {'C': float(popt[0]), 'alpha': float(popt[1])}
        print(f"\nmax_null fit: C={popt[0]:.4e}, alpha={popt[1]:.4f}")
    except Exception as e:
        print(f"max_null fit failed: {e}")

    try:
        popt2, _ = curve_fit(lambda x, C, a: C / x**a, Ls, subcancels, p0=[0.01, 1.2])
        fit_sc = {'C': float(popt2[0]), 'alpha': float(popt2[1])}
        print(f"sub_cancel fit: C={popt2[0]:.4e}, alpha={popt2[1]:.4f}")
    except Exception as e:
        print(f"sub_cancel fit failed: {e}")

    # Print dim table
    print("\n--- DIM FORMULA VERIFICATION ---")
    print(f"{'lsq':>6}  {'N':>5}  {'dim':>5}  {'exp+1':>6}  {'match':>8}  {'max_null':>10}  {'sub_cancel':>12}")
    for r in results:
        print(f"{r['lsq']:>6}  {r['N']:>5}  {r['dim']:>5}  {r['expected_plus1']:>6}  {r['formula_ok']:>8}  {r['max_null']:>10.3e}  {r['sub_cancel']:>12.3e}")

    out = {
        'thresh': THRESH, 'G': G,
        'results': results,
        'fit_maxnull': fit_mn,
        'fit_subcancel': fit_sc,
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

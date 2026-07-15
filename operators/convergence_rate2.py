#!/usr/bin/env python3
"""
convergence_rate2.py — Revised convergence rate with N = max(5*lsq, 250).

convergence_rate.py used N=3*lsq. But sc_vs_N.py showed:
- At lsq=50: dim saturates at N=4*lsq=200; sc minimum ~1.35e-4 at N=5*lsq=250
- At lsq=80: dim saturates at N=5*lsq=400; sc minimum ~5.6e-5 at saturation

This script uses N=5*lsq to ensure dim is fully saturated, then averages
over a few N values to reduce oscillation artifacts.

Tests lsq = 5, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50.
For each lsq, tests N = {4*lsq, 5*lsq, 6*lsq} and takes geometric mean.
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

OUTFILE = os.path.expanduser("~/riemann/results/convergence_rate2.json")


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


def measure_sc(lambda_sq, N_use):
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
    tau = AB - P
    tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    null_mask = np.abs(ev) < THRESH
    dim = null_mask.sum()
    if dim == 0:
        return float('nan'), int(dim), float(ev[0])
    V_null = V[:, null_mask]
    AB_proj = V_null.T @ AB @ V_null
    P_proj  = V_null.T @ P  @ V_null
    sc = float(np.linalg.norm(AB_proj - P_proj, 'fro') / (np.linalg.norm(AB_proj, 'fro') + 1e-100))
    return sc, int(dim), float(ev[0])


def main():
    lsq_list = [5, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50]
    # N multipliers: average over {4x, 5x, 6x} to reduce oscillation
    N_mults = [4, 5, 6]
    expected_dim = lambda lsq: 2 * (int(lsq) - 2) + 1

    results = []
    t0 = time.time()
    print(f"{'lsq':>5}  {'L':>5}  {'sc_avg':>12}  {'sc_min':>12}  {'dim_N5':>8}")

    for lsq in lsq_list:
        L = float(mp.log(lsq))
        exp_dim = expected_dim(lsq)
        sc_values = []
        dim_values = []
        for mult in N_mults:
            N = max(int(mult * lsq), 20)
            t1 = time.time()
            sc, dim, lam = measure_sc(lsq, N)
            dt = time.time() - t1
            sc_values.append(sc)
            dim_values.append(dim)
            print(f"  N={N} ({mult}x): sc={sc:.4e}  dim={dim}  lam={lam:.3e}  ({dt:.1f}s)", flush=True)

        sc_geomean = float(np.exp(np.mean(np.log([s for s in sc_values if not np.isnan(s)]))))
        sc_min = float(min(s for s in sc_values if not np.isnan(s)))
        dim_5x = dim_values[1] if len(dim_values) > 1 else dim_values[0]  # at 5x lsq
        print(f"  → geomean={sc_geomean:.4e}  min={sc_min:.4e}  dim={dim_5x}  (exp_dim={exp_dim})")

        results.append({
            'lsq': lsq, 'L': round(L, 6),
            'sc_values': sc_values,
            'sc_geomean': sc_geomean,
            'sc_min': sc_min,
            'dim_5x': dim_5x,
            'expected_dim': exp_dim,
            'dim_frac': dim_5x / exp_dim,
        })

    # Fit power law C/L^alpha to geomean values
    Ls = np.array([r['L'] for r in results])
    sc_gm = np.array([r['sc_geomean'] for r in results])
    sc_mn = np.array([r['sc_min'] for r in results])

    fit_gm = fit_mn = None
    try:
        popt, _ = curve_fit(lambda x, C, a: C / x**a, Ls, sc_gm, p0=[0.002, 1.5])
        fit_gm = {'C': float(popt[0]), 'alpha': float(popt[1])}
        print(f"\nGeomean fit: C={popt[0]:.4e}, alpha={popt[1]:.4f}")
    except Exception as e:
        print(f"Geomean fit failed: {e}")

    try:
        popt2, _ = curve_fit(lambda x, C, a: C / x**a, Ls, sc_mn, p0=[0.002, 1.5])
        fit_mn = {'C': float(popt2[0]), 'alpha': float(popt2[1])}
        print(f"Min fit: C={popt2[0]:.4e}, alpha={popt2[1]:.4f}")
    except Exception as e:
        print(f"Min fit failed: {e}")

    out = {
        'lsq_list': lsq_list,
        'N_mults': N_mults,
        'results': results,
        'fit_geomean': fit_gm,
        'fit_min': fit_mn,
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

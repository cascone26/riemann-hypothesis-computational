#!/usr/bin/env python3
"""
convergence_rate3.py — Extended convergence rate at lsq=60,80,100 with N=5*lsq.

Extends convergence_rate2.py to larger lsq to confirm α≈1.81 holds.
Also tests whether the oscillation structure changes at large lsq.

At N=5*lsq:
- lsq=60: N=300, matrix 601×601
- lsq=80: N=400, matrix 801×801
- lsq=100: N=500, matrix 1001×1001

Time estimates: ~5s, ~7s, ~12s per case. 3 N values each → total ~3-4 min.
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

OUTFILE = os.path.expanduser("~/riemann/results/convergence_rate3.json")


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
    # Extend to larger lsq
    lsq_list = [60, 80, 100]
    N_mults = [4, 5, 6]

    expected_dim = lambda lsq: 2 * (int(lsq) - 2) + 1

    # Also include key lsq=50 point for overlap with convergence_rate2
    lsq_list_all = [50] + lsq_list

    results = []
    t0 = time.time()

    for lsq in lsq_list_all:
        L = float(mp.log(lsq))
        exp_dim = expected_dim(lsq)
        sc_values = []
        dim_values = []
        print(f"\nlsq={lsq}, L={L:.4f}, exp_dim={exp_dim}", flush=True)
        for mult in N_mults:
            N = int(mult * lsq)
            t1 = time.time()
            sc, dim, lam = measure_sc(lsq, N)
            dt = time.time() - t1
            sc_values.append(sc)
            dim_values.append(dim)
            print(f"  N={N:4d} ({mult}x): sc={sc:.4e}  dim={dim}/{exp_dim}  lam={lam:.3e}  ({dt:.1f}s)", flush=True)

        sc_geomean = float(np.exp(np.mean(np.log([s for s in sc_values if not np.isnan(s)]))))
        sc_min = float(min(s for s in sc_values if not np.isnan(s)))
        dim_5x = dim_values[1]
        print(f"  → geomean={sc_geomean:.4e}  min={sc_min:.4e}  dim_5x={dim_5x}/{exp_dim}", flush=True)

        results.append({
            'lsq': lsq, 'L': round(L, 6),
            'sc_values': sc_values,
            'sc_geomean': sc_geomean,
            'sc_min': sc_min,
            'dim_5x': dim_5x,
            'expected_dim': exp_dim,
        })

    # Combined fit with convergence_rate2 data
    # Load convergence_rate2 results for joint fit
    cr2_file = os.path.expanduser("~/riemann/results/convergence_rate2.json")
    joint_L = []
    joint_gm = []
    try:
        with open(cr2_file) as f:
            cr2 = json.load(f)
        for r in cr2['results']:
            joint_L.append(r['L'])
            joint_gm.append(r['sc_geomean'])
        print(f"\nLoaded {len(joint_L)} points from convergence_rate2")
    except:
        print("Could not load convergence_rate2.json")

    # Add current results (lsq=60,80,100 only — lsq=50 duplicates cr2)
    for r in results:
        if r['lsq'] != 50:
            joint_L.append(r['L'])
            joint_gm.append(r['sc_geomean'])

    Ls = np.array(joint_L)
    gm = np.array(joint_gm)
    try:
        popt, pcov = curve_fit(lambda x, C, a: C / x**a, Ls, gm, p0=[0.003, 1.8])
        perr = np.sqrt(np.diag(pcov))
        print(f"\nJoint fit (lsq=5..100, N=5*lsq): C={popt[0]:.4e}±{perr[0]:.2e}, alpha={popt[1]:.4f}±{perr[1]:.4f}")
    except Exception as e:
        print(f"Joint fit failed: {e}")
        popt = [None, None]
        perr = [None, None]

    out = {
        'lsq_list': lsq_list_all,
        'N_mults': N_mults,
        'results': results,
        'joint_fit': {'C': float(popt[0]) if popt[0] else None,
                      'alpha': float(popt[1]) if popt[1] else None},
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

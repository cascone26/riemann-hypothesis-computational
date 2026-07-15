#!/usr/bin/env python3
"""
sc_vs_N.py — Does sub_cancel converge with N at fixed lsq?

Two axes of convergence:
  (a) lsq → ∞ at fixed N/lsq ratio: convergence_rate.py showed α≈1.2
  (b) N → ∞ at fixed lsq: this script tests this

If sc(0.5) decreases with N at fixed lsq → BOTH axes matter.
If sc(0.5) saturates with N → only lsq axis matters.

Tests lsq=50 at N = 50, 75, 100, 150, 200, 300, 400, 600, 800.
Also tests lsq=20 and lsq=80 for comparison.
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

OUTFILE = os.path.expanduser("~/riemann/results/sc_vs_N.json")


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
    # Primary: lsq=50 with varying N
    N_values = [50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 800]
    lsq_primary = 50

    # Also test lsq=20 and lsq=80 for comparison
    extra_cases = [
        (20, [20, 40, 60, 80, 100, 150, 200]),
        (80, [80, 120, 160, 240, 320, 400, 480]),
    ]

    results = {}
    t0 = time.time()

    print(f"=== Primary: lsq={lsq_primary} varying N ===")
    L50 = float(mp.log(lsq_primary))
    rows = []
    for N in N_values:
        ratio = N / lsq_primary
        t1 = time.time()
        sc, dim, lam_min = measure_sc(lsq_primary, N)
        dt = time.time() - t1
        rows.append({'N': N, 'ratio': round(ratio, 2), 'sc': sc, 'dim': dim, 'lam_min': lam_min, 'dt': round(dt, 2)})
        print(f"  N={N:4d} ({ratio:.1f}x): sc={sc:.4e}  dim={dim:3d}  lam={lam_min:.3e}  ({dt:.1f}s)", flush=True)
    results[str(lsq_primary)] = rows

    for lsq_extra, N_vals in extra_cases:
        print(f"\n=== lsq={lsq_extra} varying N ===")
        rows_e = []
        for N in N_vals:
            ratio = N / lsq_extra
            t1 = time.time()
            sc, dim, lam_min = measure_sc(lsq_extra, N)
            dt = time.time() - t1
            rows_e.append({'N': N, 'ratio': round(ratio, 2), 'sc': sc, 'dim': dim, 'lam_min': lam_min, 'dt': round(dt, 2)})
            print(f"  N={N:4d} ({ratio:.1f}x): sc={sc:.4e}  dim={dim:3d}  lam={lam_min:.3e}  ({dt:.1f}s)", flush=True)
        results[str(lsq_extra)] = rows_e

    out = {'results': results, 'total_time': round(time.time() - t0, 1)}
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

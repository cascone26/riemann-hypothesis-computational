#!/usr/bin/env python3
"""
dim_and_maxnull.py — Test dim formula and max_null eigenvalue scaling at large lsq.

Two studies:
1. Verify dim = 2*(floor(lsq) - 2) + 1 at lsq = 20, 50, 100, 200, 500 with N = lsq+10
2. Measure how max_null eigenvalue scales with lsq (related to sub_cancel rate)

Key question for theory: does max_null ~ 1/L^{1.232}? (matches sub_cancel rate)

Runs on Olympus. ~10 min total.
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

OUTFILE = os.path.expanduser("~/riemann/results/dim_and_maxnull.json")


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


def build_tau_eigh(lambda_sq, N_use=None):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)

    if N_use is None:
        N_use = max(64, int(lambda_sq) + 5)  # N >= lambda_sq for full formula

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
    # (lsq, N_use) pairs — N_use >= lsq for correct dim formula
    test_cases = [
        (20,  30),
        (50,  60),
        (100, 110),
        (200, 210),
        (500, 510),
    ]

    results = []
    t0 = time.time()

    for lsq, N_use in test_cases:
        print(f"\nlsq={lsq}, N={N_use}...", flush=True)
        t1 = time.time()
        ev, V, idx, L, AB, P = build_tau_eigh(lsq, N_use)
        dt = time.time() - t1

        null_mask = np.abs(ev) < THRESH
        dim = null_mask.sum()
        lam_min = ev[0]
        max_null = ev[null_mask].max() if dim > 0 else float('nan')

        # Expected dim formula
        expected_dim = 2 * (int(lsq) - 2) + 1

        # Sub_cancel
        V_null = V[:, null_mask]
        AB_proj = V_null.T @ AB @ V_null
        P_proj  = V_null.T @ P  @ V_null
        sub_cancel = np.linalg.norm(AB_proj - P_proj, 'fro') / (np.linalg.norm(AB_proj, 'fro') + 1e-100)

        print(f"  dim={dim} (expected {expected_dim}, {'OK' if dim == expected_dim else 'MISMATCH'})")
        print(f"  lam_min={lam_min:.4e}, max_null={max_null:.4e}")
        print(f"  sub_cancel={sub_cancel:.4e}  L={L:.4f}  ({dt:.1f}s)")

        results.append({
            'lsq': lsq, 'N': N_use, 'L': round(L, 6),
            'dim': int(dim), 'expected_dim': expected_dim,
            'formula_ok': int(dim == expected_dim),
            'lam_min': round(float(lam_min), 12),
            'max_null': round(float(max_null), 12),
            'sub_cancel': round(float(sub_cancel), 12),
        })

    # Fit power law for max_null ~ C / L^alpha
    from scipy.optimize import curve_fit
    Ls = np.array([r['L'] for r in results])
    maxnull = np.array([r['max_null'] for r in results])
    subcancel = np.array([r['sub_cancel'] for r in results])

    print("\n--- Power law fit: max_null ~ C / L^alpha ---")
    try:
        popt, _ = curve_fit(lambda x, C, a: C / x**a, Ls, maxnull, p0=[0.01, 1.2])
        C_mn, a_mn = popt
        print(f"  max_null: C={C_mn:.4e}, alpha={a_mn:.4f}")
    except Exception as e:
        print(f"  Fit failed: {e}")
        C_mn, a_mn = None, None

    print("\n--- Power law fit: sub_cancel ~ C / L^alpha ---")
    try:
        popt2, _ = curve_fit(lambda x, C, a: C / x**a, Ls, subcancel, p0=[0.01, 1.2])
        C_sc, a_sc = popt2
        print(f"  sub_cancel: C={C_sc:.4e}, alpha={a_sc:.4f}")
    except Exception as e:
        print(f"  Fit failed: {e}")
        C_sc, a_sc = None, None

    out = {
        'thresh': THRESH, 'G': G,
        'results': results,
        'fit_maxnull': {'C': C_mn, 'alpha': a_mn} if C_mn else None,
        'fit_subcancel': {'C': C_sc, 'alpha': a_sc} if C_sc else None,
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

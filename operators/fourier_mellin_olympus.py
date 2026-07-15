#!/usr/bin/env python3
"""
fourier_mellin_olympus.py — Large-scale Fourier-Mellin zero extraction on Olympus.

Runs at multiple (lambda_sq, N) pairs to see how zero detection improves.
Uses min_kinetic eigenvector selection. Reports |F(gamma_k)|/max at all known zeros.

Designed for Olympus (fast CPU, no GPU needed for this task).
Output: ~/riemann/results/fourier_mellin_olympus.json
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp
from scipy.signal import argrelmin

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000

KNOWN_ZEROS = [
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302,
]

OUTFILE = os.path.expanduser("~/riemann/results/fourier_mellin_olympus.json")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)


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
    return (tau + tau.T) / 2.0, idx.astype(int)


def fourier_mellin(c, idx, L, z_grid):
    alpha = TWO_PI * idx[:, None] - L * z_grid[None, :]
    half  = 0.5 * alpha
    with np.errstate(divide='ignore', invalid='ignore'):
        kern = np.where(np.abs(half) < 1e-10, 1.0, np.sin(half) / half)
    return np.sqrt(L) * np.sum(c[:, None] * np.exp(1j * half) * kern, axis=0)


def min_kinetic_eigvec(ev, evec, idx, thresh=1e-2):
    mask = np.abs(ev) < thresh
    nn_evec = evec[:, mask]
    K = int(mask.sum())
    if K == 0:
        return evec[:, 0], 0
    n2 = (idx ** 2).astype(float)
    K_proj = (nn_evec * n2[:, None]).T @ nn_evec
    K_proj = (K_proj + K_proj.T) / 2
    kev, kevec = np.linalg.eigh(K_proj)
    c = nn_evec @ kevec[:, 0]
    return c / np.linalg.norm(c), K


def run_case(lambda_sq, N, z_max=60.0, Nz=200000, label=""):
    L = float(np.log(lambda_sq))
    print(f"\n{'='*65}")
    print(f"lambda_sq={lambda_sq}, N={N}, L={L:.4f} {label}", flush=True)

    t0 = time.time()
    tau, idx = build_tau(lambda_sq, N)
    t1 = time.time()
    ev, evec = scipy.linalg.eigh(tau)
    t2 = time.time()
    print(f"  build={t1-t0:.1f}s  eigh={t2-t1:.1f}s  lam_min={ev[0]:.4e}", flush=True)

    K_02 = int((np.abs(ev) < 1e-2).sum())
    K_03 = int((np.abs(ev) < 1e-3).sum())
    print(f"  K(1e-2)={K_02}  K(1e-3)={K_03}", flush=True)

    c_kin, K = min_kinetic_eigvec(ev, evec, idx, thresh=1e-2)

    z_grid = np.linspace(5.0, z_max, Nz)
    F      = fourier_mellin(c_kin, idx, L, z_grid)
    F_abs  = np.abs(F)
    Fmax   = F_abs.max()
    t3 = time.time()
    print(f"  F-M computed in {t3-t2:.1f}s", flush=True)

    # Check at known zeros
    zero_checks = []
    print(f"  {'k':>3}  {'gamma_k':>10}  {'|F|/max':>12}  {'local_min':>10}  {'detected_z':>12}  {'error':>10}")
    for ki, gk in enumerate(KNOWN_ZEROS):
        iz = np.argmin(np.abs(z_grid - gk))
        fval = float(F_abs[iz])
        lo, hi = max(0, iz-300), min(len(F_abs)-1, iz+300)
        is_min = bool(fval == F_abs[lo:hi+1].min())

        # Local refinement: find the minimum in a window around gamma_k
        window = 0.3
        wlo = np.searchsorted(z_grid, gk - window)
        whi = np.searchsorted(z_grid, gk + window)
        if whi > wlo:
            local_imin = wlo + np.argmin(F_abs[wlo:whi])
            z_detected = float(z_grid[local_imin])
            f_detected = float(F_abs[local_imin])
        else:
            z_detected = float(z_grid[iz])
            f_detected = fval

        error = z_detected - gk
        mark = '*' if is_min else ' '
        print(f"  {ki+1:>3}  {gk:>10.5f}  {fval/Fmax:>12.4e}  {mark:>10}  {z_detected:>12.6f}  {error:>+10.6f}", flush=True)

        zero_checks.append({
            'k': ki + 1,
            'gamma_true': gk,
            'F_abs': fval,
            'F_abs_rel': fval / Fmax,
            'is_local_min': is_min,
            'z_detected': z_detected,
            'z_detected_Fabs_rel': f_detected / Fmax,
            'error': error,
        })

    return {
        'lambda_sq': lambda_sq, 'N': N, 'L': L, 'label': label,
        'K': K_02, 'lam_min': float(ev[0]),
        'zero_checks': zero_checks,
        'dt': round(t3 - t0, 1),
    }


def main():
    cases = [
        # (lambda_sq, N, label)
        (50, 192, "baseline"),
        (50, 400, "high_N"),
        (100, 300, "lambda100"),
        (200, 400, "lambda200"),
        (500, 500, "lambda500"),
    ]

    all_results = []
    for lsq, N, label in cases:
        try:
            r = run_case(lsq, N, label=label)
            all_results.append(r)
            with open(OUTFILE, 'w') as f:
                json.dump(all_results, f, indent=2)
            print(f"  -> saved {len(all_results)}/{len(cases)}", flush=True)
        except Exception as e:
            print(f"ERROR lsq={lsq} N={N}: {e}", flush=True)
            import traceback; traceback.print_exc()


if __name__ == "__main__":
    main()

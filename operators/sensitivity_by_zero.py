#!/usr/bin/env python3
"""
sensitivity_by_zero.py — Sensitivity |dlam/dσ| for each of the first 15 Riemann zeros.

At lsq=50 (N=250), measure how much EACH zero contributes to sensitivity.
Higher zeros (larger γ) are further from the prime range → smaller sensitivity?
Or does each zero contribute ~1.0 per unit Δσ regardless of γ?

Also tests lsq=250 for γ₁..γ₅ (larger range → more zeros enter the prime range).
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
DELTA_SIGMA = 0.002
OUTFILE = os.path.expanduser("~/riemann/results/sensitivity_by_zero.json")

RIEMANN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831778, 65.112544,
]


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


def build_tau_components(lsq, N_use):
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lsq)
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
    return arch, B, P, E, grid, w, L, ns


def zero_delta_A(E, grid, w, sigma_new, gamma, sigma_old=0.5):
    s0 = 0.25 + 0.5j * grid
    def contrib(sigma):
        rho = sigma + 1j * gamma
        rho_bar = sigma - 1j * gamma
        return (np.real(1.0 / (s0 - rho) + 1.0 / rho) +
                np.real(1.0 / (s0 - rho_bar) + 1.0 / rho_bar))
    delta_Om = contrib(sigma_new) - contrib(sigma_old)
    return ((E * (w * delta_Om)) @ E.conj().T).real / (2.0 * np.pi)


def measure_sensitivity(arch, B, P, E, grid, w, gamma, lam0, ds=DELTA_SIGMA):
    dA = zero_delta_A(E, grid, w, 0.5 + ds, gamma)
    tau2 = arch + B + dA - P
    tau2 = (tau2 + tau2.T) / 2.0
    lam_g = float(scipy.linalg.eigvalsh(tau2, subset_by_index=[0, 0])[0])
    slope = (lam_g - lam0) / ds
    return lam_g, float(slope)


def main():
    results_all = {}
    t0 = time.time()

    # Test at lsq=50
    for lsq, nzeros in [(50, 15), (250, 8)]:
        N = 5 * lsq
        L = float(mp.log(lsq))
        print(f"\n=== lsq={lsq}, N={N} ===", flush=True)
        arch, B, P, E, grid, w, L, ns = build_tau_components(lsq, N)
        tau = arch + B - P
        tau = (tau + tau.T) / 2.0
        lam0 = float(scipy.linalg.eigvalsh(tau, subset_by_index=[0, 0])[0])
        print(f"lam0 = {lam0:.4e}", flush=True)

        results = []
        for k, gamma_k in enumerate(RIEMANN_ZEROS[:nzeros]):
            t1 = time.time()
            lam_g, slope = measure_sensitivity(arch, B, P, E, grid, w, gamma_k, lam0)
            dt = time.time() - t1
            inside = "in" if float(np.sqrt(lsq)) > gamma_k else "out"
            print(f"  k={k+1:2d} γ={gamma_k:8.3f} ({inside}): slope={slope:.4e}  lam_ghost={lam_g:.4e}  ({dt:.1f}s)", flush=True)
            results.append({
                'k': k + 1, 'gamma': gamma_k,
                'inside_prime_range': float(np.sqrt(lsq)) > gamma_k,
                'lam_ghost': lam_g,
                'dlam_dsigma': float(slope),
                'dt': round(dt, 1),
            })
        results_all[str(lsq)] = {'lam0': lam0, 'L': round(L, 4), 'results': results}

    out = {
        'delta_sigma': DELTA_SIGMA,
        'riemann_zeros_used': RIEMANN_ZEROS[:15],
        'data': results_all,
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
sensitivity_large.py — Sensitivity scaling at lsq=100,150,200,250.

Key threshold: γ₁=14.135 enters the prime range when λ=√lsq≈14.135 → lsq≈200.
At lsq≥200, the zero γ₁ is "inside" the prime sum. Does sensitivity grow further?

Tests lsq=100,150,200,250 to extend sensitivity_scaling.py.
N=4*lsq (slightly below 5* to save time at large lsq; close enough to saturation).
Only Δσ=0.001 for speed (confirmed linear).
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
GAMMA1 = 14.134725
DELTA_SIGMA = 0.001
OUTFILE = os.path.expanduser("~/riemann/results/sensitivity_large.json")


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


def build_and_measure(lsq, N_use, ds):
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
    AB = arch + B
    # Baseline
    tau = (AB - P); tau = (tau + tau.T) / 2.0
    lam0 = float(scipy.linalg.eigvalsh(tau, subset_by_index=[0, 0])[0])
    # Sub-cancel
    ev, V = scipy.linalg.eigh(tau)
    null_mask = np.abs(ev) < THRESH
    dim = null_mask.sum()
    if dim > 0:
        Vn = V[:, null_mask]
        ABp = Vn.T @ AB @ Vn
        Pp  = Vn.T @ P  @ Vn
        sc = float(np.linalg.norm(ABp - Pp, 'fro') / (np.linalg.norm(ABp, 'fro') + 1e-100))
    else:
        sc = float('nan')
    # Moved zero
    s0 = 0.25 + 0.5j * grid
    def contrib(sigma):
        rho = sigma + 1j * GAMMA1
        rho_bar = sigma - 1j * GAMMA1
        return (np.real(1.0 / (s0 - rho) + 1.0 / rho) +
                np.real(1.0 / (s0 - rho_bar) + 1.0 / rho_bar))
    dOm = contrib(0.5 + ds) - contrib(0.5)
    dA = ((E * (w * dOm)) @ E.conj().T).real / (2.0 * np.pi)
    AB2 = AB + dA
    tau2 = (AB2 - P); tau2 = (tau2 + tau2.T) / 2.0
    lam_ghost = float(scipy.linalg.eigvalsh(tau2, subset_by_index=[0, 0])[0])
    slope = (lam_ghost - lam0) / ds
    return lam0, sc, int(dim), lam_ghost, float(slope)


def main():
    lsq_list = [100, 150, 200, 250]
    results = []
    t0 = time.time()

    for lsq in lsq_list:
        N = 4 * lsq  # slightly below optimal 5× but saves time
        L = float(mp.log(lsq))
        lambda_val = float(np.sqrt(lsq))
        inside = "INSIDE" if lambda_val > GAMMA1 else "outside"
        print(f"\nlsq={lsq}, N={N}, L={L:.4f}, λ={lambda_val:.2f} ({inside} γ₁=14.135)", flush=True)
        t1 = time.time()
        lam0, sc, dim, lam_ghost, slope = build_and_measure(lsq, N, DELTA_SIGMA)
        dt = time.time() - t1
        print(f"  lam0={lam0:.4e}  sc={sc:.4e}  dim={dim}", flush=True)
        print(f"  lam_ghost(Δσ=0.001)={lam_ghost:.4e}  slope={slope:.4e}  ({dt:.1f}s)", flush=True)
        results.append({
            'lsq': lsq, 'L': round(L, 6), 'N': N,
            'lambda': round(lambda_val, 4),
            'gamma1_inside': lambda_val > GAMMA1,
            'lam0': lam0, 'sc': sc, 'dim': dim,
            'lam_ghost': lam_ghost,
            'dlam_dsigma': slope,
            'dt': round(dt, 1),
        })

    out = {'lsq_list': lsq_list, 'gamma1': GAMMA1, 'delta_sigma': DELTA_SIGMA,
           'results': results, 'total_time': round(time.time() - t0, 1)}
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

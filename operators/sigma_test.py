#!/usr/bin/env python3
"""
sigma_test.py — Test sub_cancel convergence for σ ≠ 1/2 in prime weight.

The key claim: sub_cancel → 0 as λ²→∞ ONLY when σ=1/2 (on critical line).

We test: replace Λ(k)/k^{1/2} with Λ(k)/k^σ for σ = 0.3, 0.4, 0.5, 0.6, 0.7.
If the claim is correct:
- σ=0.5: sub_cancel → 0 (standard RH case)
- σ≠0.5: sub_cancel → constant > 0 (or grows)

This is the NUMERICAL FALSIFICATION TEST for "sub_cancel → 0 ⟺ σ=1/2".

Note: this doesn't directly prove/disprove RH, but shows that the cancellation
is SENSITIVE to the critical exponent σ=1/2 — a necessary condition for our
approach to be connected to RH.

Runs in ~5 min on Olympus. Uses N=3*lsq for proper sub_cancel measurement.
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

OUTFILE = os.path.expanduser("~/riemann/results/sigma_test.json")


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


def build_tau_sigma(lambda_sq, N_use, sigma):
    """Build τ_σ = A + B - P_σ where P_σ uses Λ(k)/k^σ instead of Λ(k)/k^{1/2}."""
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
        # Key change: use k^(-sigma) instead of k^(-1/2)
        # VM has (k, log_p). We need k^(-sigma).
        # Standard: wts = Λ(k) × k^(-1/2) = log(prime_base) × k^(-1/2)
        # Modified: wts = Λ(k) × k^(-sigma) = log(prime_base) × k^(-sigma)
        vm_k = np.array([k for k, _ in VM])
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
    V_null = V[:, null_mask]
    if dim == 0:
        return {'dim': 0, 'lam_min': float(ev[0]), 'sub_cancel': float('nan')}
    AB_proj = V_null.T @ AB @ V_null
    P_proj  = V_null.T @ P  @ V_null
    sc = np.linalg.norm(AB_proj - P_proj, 'fro') / (np.linalg.norm(AB_proj, 'fro') + 1e-100)
    return {
        'dim': int(dim),
        'lam_min': float(ev[0]),
        'sub_cancel': float(sc),
        'max_null': float(ev[null_mask].max()) if dim > 0 else float('nan'),
    }


def main():
    # Test at lsq = 5, 10, 20, 30, 50 for sigma = 0.3, 0.4, 0.5, 0.6, 0.7
    lsq_list = [5, 10, 20, 30, 50, 80]
    sigma_list = [0.3, 0.4, 0.5, 0.6, 0.7]
    N_factor = 3  # N = 3 * lsq

    all_results = {}
    t0 = time.time()

    for sigma in sigma_list:
        print(f"\n=== sigma={sigma} ===", flush=True)
        results = []
        for lsq in lsq_list:
            N_use = int(np.ceil(N_factor * lsq))
            L = float(mp.log(lsq))
            t1 = time.time()
            r = build_tau_sigma(lsq, N_use, sigma)
            dt = time.time() - t1
            print(f"  lsq={lsq:3d} N={N_use:3d}: dim={r['dim']:3d}  sc={r['sub_cancel']:.4e}  lam_min={r['lam_min']:.3e}  ({dt:.1f}s)", flush=True)
            r.update({'lsq': lsq, 'N': N_use, 'L': round(L, 4)})
            results.append(r)
        all_results[str(sigma)] = results

    # Print comparison table
    print("\n--- sub_cancel by sigma and lsq ---")
    header = f"{'lsq':>5}  {'L':>5}" + "".join(f"  σ={s:.1f}" for s in sigma_list)
    print(header)
    for i, lsq in enumerate(lsq_list):
        L = float(mp.log(lsq))
        row = f"{lsq:>5}  {L:>5.3f}"
        for sigma in sigma_list:
            sc = all_results[str(sigma)][i]['sub_cancel']
            row += f"  {sc:.3e}"
        print(row)

    out = {
        'lsq_list': lsq_list,
        'sigma_list': sigma_list,
        'N_factor': N_factor,
        'thresh': THRESH,
        'results': all_results,
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

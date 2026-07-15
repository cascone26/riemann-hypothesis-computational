#!/usr/bin/env python3
"""
eigvec_convergence.py — Track how the minimum eigenvector of tau_N evolves with N.

Key question: as N grows, does the minimum eigenvector converge to a specific
function? If so, which Fourier components (n values) dominate?

The limit function is the object we need to prove positivity for analytically.
If it's purely the n=0 mode (constant), that's a tractable target.
"""

import json, os
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30

LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/eigvec_convergence.json")


def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
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

VM = vm_table(LAMBDA_SQ)


def build_tau(N):
    idx = list(range(-N, N+1))
    D = len(idx)
    R = TWO_PI * N / L + 80.0
    G = max(60_000, 3000 * N)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5*h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI

    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        s = TWO_PI * n / L
        d = s - grid
        sm = np.abs(d) < 1e-12
        E[k, sm] = np.sqrt(L)
        dd = d[~sm]
        E[k, ~sm] = (1/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2*np.sin(0.5*dd*L) / dd

    arch = ((E * (w * Om)) @ E.conj().T).real / (2*np.pi)

    sh = np.sinh(L/4)**2
    ns = np.array(idx, dtype=float)
    B = 32*L*sh*(L**2 - 16*np.pi**2*np.outer(ns,ns)) / np.outer(L**2+16*np.pi**2*ns**2, L**2+16*np.pi**2*ns**2)

    log_k = np.array([np.log(k) for k, _ in VM])
    wts = np.array([lk * k**(-0.5) for k, lk in VM])
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            q = np.where(n==m,
                         2*(1-np.abs(log_k)/L)*np.cos(TWO_PI*n*log_k/L),
                         (np.sin(TWO_PI*m*log_k/L)-np.sin(TWO_PI*n*log_k/L))/(np.pi*(n-m)+1e-300))
            P[ii,jj] = P[jj,ii] = float(wts @ q)

    tau = arch + B - P
    return (tau + tau.T)/2, idx


def analyze_eigvec(N):
    print(f"N={N}...", flush=True)
    tau, idx = build_tau(N)
    ev, evec = np.linalg.eigh(tau)
    v_min = evec[:, 0]   # minimum eigenvector in idx order

    # Power in each n-mode: |v[n]|^2, normalized
    power = np.abs(v_min)**2
    power /= power.sum()

    # Which n-modes carry most power?
    top_k = 10
    top_idx = np.argsort(power)[::-1][:top_k]
    top_modes = [(idx[i], float(power[i])) for i in top_idx]

    # Power in n=0 mode
    n0_idx = idx.index(0)
    p0 = float(power[n0_idx])

    # Power in |n|≤5 vs rest
    low_n_mask = np.array([abs(n) <= 5 for n in idx])
    p_low = float(power[low_n_mask].sum())

    # Cumulative power in |n|≤k for k=0..10
    cum_power = []
    for k in range(11):
        mask = np.array([abs(n) <= k for n in idx])
        cum_power.append(float(power[mask].sum()))

    print(f"  lam_min={ev[0]:.4e}  p(n=0)={p0:.4f}  p(|n|≤5)={p_low:.4f}", flush=True)
    print(f"  Top 5 modes: {top_modes[:5]}", flush=True)

    return {
        "N": N,
        "lambda_min": float(ev[0]),
        "power_n0": p0,
        "power_abs_n_le5": p_low,
        "cumulative_power_abs_n_le_k": cum_power,  # k=0..10
        "top10_modes": top_modes,
        "5_smallest_ev": ev[:5].tolist()
    }


def main():
    print("EIGENVECTOR CONVERGENCE ANALYSIS", flush=True)
    print(f"lambda^2={LAMBDA_SQ}, L={L:.4f}", flush=True)

    N_values = [16, 24, 32, 48, 64, 96, 128]
    results = []

    for N in N_values:
        r = analyze_eigvec(N)
        results.append(r)
        with open(OUTFILE, "w") as f:
            json.dump({"lambda_sq": LAMBDA_SQ, "L": L, "results": results}, f, indent=2)

    print("\n=== CONVERGENCE SUMMARY ===", flush=True)
    print(f"{'N':>4}  {'lam_min':>10}  {'p(n=0)':>8}  {'p(|n|≤1)':>10}  {'p(|n|≤5)':>10}", flush=True)
    for r in results:
        cp = r["cumulative_power_abs_n_le_k"]
        print(f"  {r['N']:>4}  {r['lambda_min']:>10.3e}  {r['power_n0']:>8.4f}  "
              f"{cp[1]:>10.4f}  {cp[5]:>10.4f}", flush=True)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()

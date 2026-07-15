#!/usr/bin/env python3
"""
tau_subspace_analysis.py — Analyze the near-zero subspace of tau_N.

Key questions:
1. How many eigenvalues of tau_N are within epsilon of zero?
   (Should match count of zeta zeros with gamma < bandwidth s_N = 2*pi*N/L)
2. What is tau[0,0] explicitly? If > 0, n=0 mode is NOT in the zero subspace.
3. What is the quadratic form for each basis mode e_n separately?
   (Diagonal of tau: tau[n,n] for each n)
4. Does the zero subspace grow as N grows, and at what rate?

The diagonal tau[n,n] = W(e_n, e_n) tells us which directions are "hard"
(near zero) and which are easily positive.
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
                       "../results/tau_subspace_analysis.json")

# Reference zeros (first 50 if computable)
print("Computing reference zeros...", flush=True)
REFZEROS = [float(mp.im(mp.zetazero(k))) for k in range(1, 51)]
print(f"  gamma_1..5: {[round(z,3) for z in REFZEROS[:5]]}", flush=True)


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
            q = np.where(n == m,
                         2*(1-np.abs(log_k)/L)*np.cos(TWO_PI*n*log_k/L),
                         (np.sin(TWO_PI*m*log_k/L)-np.sin(TWO_PI*n*log_k/L))/(np.pi*(n-m)+1e-300))
            P[ii,jj] = P[jj,ii] = float(wts @ q)

    tau = arch + B - P
    return (tau + tau.T)/2, idx


def analyze_N(N):
    print(f"\n=== N={N} ===", flush=True)
    tau, idx = build_tau(N)
    ev = np.linalg.eigvalsh(tau)

    # Bandwidth ceiling
    s_N = TWO_PI * N / L
    zeros_in_bw = [z for z in REFZEROS if z < s_N]
    n_zeros_bw = len(zeros_in_bw)

    # Near-zero eigenvalue count at various thresholds
    thresholds = [1e-2, 1e-3, 1e-4, 1e-5, 1e-8, 1e-10]
    near_zero = {f"{t:.0e}": int(np.sum(np.abs(ev) < t)) for t in thresholds}

    # Diagonal of tau: tau[n,n] for each n in idx
    diag = tau.diagonal()
    # The n=0 element
    n0_idx = idx.index(0)
    tau_00 = float(diag[n0_idx])

    # All diagonal elements sorted by |value|
    diag_sorted_by_abs = sorted(zip(idx, diag.tolist()), key=lambda x: abs(x[1]))

    print(f"  s_N = 2pi*N/L = {s_N:.3f}", flush=True)
    print(f"  zeros with gamma < s_N: {n_zeros_bw} (zeros: {[round(z,2) for z in zeros_in_bw[:5]]}...)", flush=True)
    print(f"  near-zero counts: {near_zero}", flush=True)
    print(f"  tau[0,0] = {tau_00:.8e}", flush=True)
    print(f"  lambda_min = {ev[0]:.6e}, lambda_max = {ev[-1]:.6e}", flush=True)
    print(f"  5 smallest ev: {ev[:5].tolist()}", flush=True)
    print(f"  10 smallest |diag| (n, tau[n,n]):", flush=True)
    for n, v in diag_sorted_by_abs[:10]:
        print(f"    n={n:>4}: tau[n,n] = {v:.6e}", flush=True)

    return {
        "N": N,
        "s_N": float(s_N),
        "zeros_in_bandwidth": n_zeros_bw,
        "lambda_min": float(ev[0]),
        "lambda_max": float(ev[-1]),
        "5_smallest_ev": ev[:5].tolist(),
        "near_zero_counts": near_zero,
        "tau_00": tau_00,
        "tau_diagonal_all": [(int(n), float(v)) for n, v in zip(idx, diag.tolist())],
        "10_smallest_abs_diag": [(int(n), float(v)) for n, v in diag_sorted_by_abs[:10]],
    }


def main():
    print("TAU SUBSPACE ANALYSIS", flush=True)
    print(f"lambda^2={LAMBDA_SQ}, L={L:.4f}", flush=True)

    N_values = [16, 32, 48, 64, 96, 128]
    results = []

    for N in N_values:
        r = analyze_N(N)
        results.append(r)
        with open(OUTFILE, "w") as f:
            json.dump({"lambda_sq": LAMBDA_SQ, "L": L, "results": results}, f, indent=2)

    print("\n=== CROSS-N SUMMARY ===", flush=True)
    print(f"{'N':>4}  {'s_N':>7}  {'zeros_in_bw':>11}  {'near0_1e-4':>10}  {'near0_1e-10':>12}  {'tau[0,0]':>12}", flush=True)
    for r in results:
        print(f"  {r['N']:>4}  {r['s_N']:>7.2f}  {r['zeros_in_bandwidth']:>11}  "
              f"{r['near_zero_counts']['1e-04']:>10}  {r['near_zero_counts']['1e-10']:>12}  "
              f"{r['tau_00']:>12.4e}", flush=True)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()

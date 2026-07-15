#!/usr/bin/env python3
"""
limit_function_extrap.py — Characterize the limit of the minimum eigenvector.

The minimum eigenvector g_min(N) of tau_N achieves lambda_min ~ machine precision.
As N->inf, g_min converges to a limit function g_lim.

Questions:
  1. Does the coefficient vector c(N) = g_min in mode basis converge? (Cauchy test)
  2. What is the limit function's mode profile? (which n modes dominate)
  3. Richardson-extrapolate lambda_min(N): does it -> 0 exactly or a positive const?
  4. Compute B - |A| - P for the limit: is the cancellation EXACT (a theorem) or approximate?

Key: if tau_inf has a genuine null vector (lambda=0 exactly), the operator is PSD
with kernel — the RH-consistent boundary case. If lambda_min -> c > 0, that would
be even stronger (strict positivity). If lambda_min -> c < 0, RH would be FALSE
(but all evidence says c = 0 from below at machine precision).
"""

import json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
LAMBDA_SQ = 14
L = float(mp.log(LAMBDA_SQ))
sh = float(np.sinh(L/4)**2)

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/limit_function_extrap.json")


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


def build_tau(N):
    VM = vm_table(LAMBDA_SQ)
    idx = list(range(-N, N+1))
    D = len(idx)
    R = TWO_PI * N / L + 80.0
    G = 60000
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
    return (tau + tau.T)/2, arch, B, P, idx


def main():
    Ns = [24, 32, 48, 64, 96, 128]
    data = []
    prev_coeffs = None

    print("Limit function extrapolation (min eigenvector as N grows)")
    print(f"{'N':>4}  {'lambda_min':>13}  {'2nd_smallest':>13}  {'coeff_shift':>12}  {'A':>8}  {'B':>8}  {'P':>8}  {'B-|A|-P':>12}")

    for N in Ns:
        tau, arch, B, P, idx = build_tau(N)
        ev, evec = np.linalg.eigh(tau)
        g = evec[:, 0]  # minimum eigenvector
        lam = ev[0]

        # A, B, P for this eigenvector
        A_val = float(g @ arch @ g)
        B_val = float(g @ B @ g)
        P_val = float(g @ P @ g)
        cancellation = B_val - abs(A_val) - P_val  # should be ~ lambda_min if A<0

        # Coefficient profile centered on n=0: store c[n] for |n|<=20
        center = N
        core = {}
        for n in range(-20, 21):
            i = center + n
            if 0 <= i < len(g):
                core[n] = float(g[i])
        # sign-normalize (make c[0] positive)
        if core.get(0, 0) < 0:
            core = {k: -v for k, v in core.items()}
            g = -g

        # Compare to previous N (only common core modes)
        coeff_shift = None
        if prev_coeffs is not None:
            common = [n for n in core if n in prev_coeffs]
            diff = np.array([core[n] - prev_coeffs[n] for n in common])
            coeff_shift = float(np.linalg.norm(diff))
        prev_coeffs = core

        # 2nd smallest for gap reference
        ev2 = ev[1]

        data.append({
            'N': N, 'lambda_min': float(lam), 'lambda_2': float(ev2),
            'A': A_val, 'B': B_val, 'P': P_val, 'cancellation': cancellation,
            'core_coeffs': core, 'coeff_shift': coeff_shift
        })

        cs = f"{coeff_shift:.4e}" if coeff_shift is not None else "  --"
        print(f"{N:>4}  {lam:>13.5e}  {ev2:>13.5e}  {cs:>12}  {A_val:>8.4f}  {B_val:>8.4f}  {P_val:>8.4f}  {cancellation:>12.5e}")

    # Richardson extrapolation of A, B, P values (they converge, unlike lambda which is ~0)
    print("\n=== Convergence of A, B, P (the limit function values) ===")
    for key in ['A', 'B', 'P']:
        vals = [d[key] for d in data]
        # fit v(N) = v_inf + c/N^p
        Nsarr = np.array([d['N'] for d in data], dtype=float)
        v = np.array(vals)
        # estimate v_inf via 3-point Richardson using last 3
        v1, v2, v3 = v[-3], v[-2], v[-1]
        N1, N2, N3 = Nsarr[-3], Nsarr[-2], Nsarr[-1]
        # assume v = v_inf + c*N^-1
        # v_inf ~ v3 + (v3-v2)*N3/(N2-N3) ... use simple extrapolation
        # Aitken-style
        try:
            denom = (v3 - v2) - (v2 - v1)
            v_inf = v3 - (v3 - v2)**2 / denom if abs(denom) > 1e-15 else v3
        except Exception:
            v_inf = v3
        print(f"  {key}_inf ≈ {v_inf:.6f}  (last value {v3:.6f})")

    # The key extrapolation: does |A|+P -> B exactly?
    print("\n=== Cancellation identity: B_inf vs |A_inf|+P_inf ===")
    A_inf = data[-1]['A']; B_inf = data[-1]['B']; P_inf = data[-1]['P']
    print(f"  B_inf = {B_inf:.6f}")
    print(f"  |A_inf| + P_inf = {abs(A_inf) + P_inf:.6f}")
    print(f"  Difference (lambda_min proxy) = {B_inf - abs(A_inf) - P_inf:.6e}")

    # Coefficient convergence -> is the limit function well-defined?
    print("\n=== Coefficient vector convergence (Cauchy test) ===")
    shifts = [d['coeff_shift'] for d in data if d['coeff_shift'] is not None]
    print(f"  Coeff shifts between consecutive N: {[f'{s:.4e}' for s in shifts]}")
    if len(shifts) >= 2:
        print(f"  Shrinking? {shifts[-1] < shifts[0]}  (ratio {shifts[-1]/shifts[0]:.3f})")

    # Print the limit function's mode profile
    print("\n=== Limit function mode profile (N=128, |n|<=12) ===")
    core = data[-1]['core_coeffs']
    for n in range(-12, 13):
        if n in core:
            bar = '#' * int(abs(core[n]) * 50)
            print(f"  n={n:>3}: {core[n]:>9.5f}  {bar}")

    with open(OUTFILE, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"\nResults: {OUTFILE}")


if __name__ == "__main__":
    main()

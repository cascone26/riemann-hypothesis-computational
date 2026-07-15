#!/usr/bin/env python3
"""
spectral_gap_vs_lambda.py — How does the spectral gap depend on lambda^2?

For each lambda^2, compute at N=32:
  - near-zero count (threshold 1e-4)
  - first eigenvalue above 1e-4 (the spectral gap)
  - tau[0,0]
  - lambda_min

Does the gap scale with number of prime powers? With lambda^2?
"""

import json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/spectral_gap_vs_lambda.json")


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


def build_tau(lambda_sq, N=32):
    L = float(mp.log(lambda_sq))
    sh = float(np.sinh(L/4)**2)
    VM = vm_table(lambda_sq)
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
    return (tau + tau.T)/2, len(VM)


def main():
    N = 32
    # Sample a range of lambda^2 values
    lambda_sq_values = [
        4, 6, 9, 14, 16, 20, 25, 30, 36, 49, 50,
        64, 81, 100, 121, 144, 169, 196, 225, 250, 289, 361, 400
    ]

    results = []
    print(f"Spectral gap vs lambda^2 (N={N})", flush=True)
    print(f"{'lsq':>6}  {'#primes':>7}  {'tau00':>8}  {'lmin':>12}  {'nz':>4}  {'gap':>12}  {'t':>5}", flush=True)

    for lsq in lambda_sq_values:
        t0 = time.time()
        try:
            tau, n_pp = build_tau(lsq, N)
            ev = np.linalg.eigvalsh(tau)
            t1 = time.time()

            nz = int(np.sum(np.abs(ev) < 1e-4))
            above = ev[ev > 1e-4]
            gap = float(above[0]) if len(above) > 0 else None

            L = float(np.log(lsq))
            r = {
                'lambda_sq': lsq, 'L': L, 'n_prime_powers': n_pp,
                'tau_00': float(tau[N, N]),
                'lambda_min': float(ev[0]), 'lambda_max': float(ev[-1]),
                'near_zero_count': nz,
                'first_gap': gap,
                'smallest_10': ev[:10].tolist(),
                'time_s': round(t1-t0, 2)
            }
            results.append(r)

            gap_str = f"{gap:.4e}" if gap else "None"
            print(f"  {lsq:>5}  {n_pp:>7}  {tau[N,N]:>8.5f}  {ev[0]:>12.4e}  "
                  f"{nz:>4}  {gap_str:>12}  {t1-t0:>5.2f}s", flush=True)

        except Exception as e:
            print(f"  {lsq:>5}: ERROR {e}", flush=True)
            continue

        with open(OUTFILE, 'w') as f:
            json.dump(results, f, indent=2)

    # Summary
    print("\n=== Summary ===", flush=True)
    print(f"{'lsq':>6}  {'n_pp':>4}  {'nz':>4}  {'gap':>12}  {'gap/lsq':>12}", flush=True)
    for r in results:
        gap_str = f"{r['first_gap']:.4e}" if r['first_gap'] else "None"
        gap_norm = f"{r['first_gap']/r['lambda_sq']:.4e}" if r['first_gap'] else "N/A"
        print(f"  {r['lambda_sq']:>5}  {r['n_prime_powers']:>4}  {r['near_zero_count']:>4}  {gap_str:>12}  {gap_norm:>12}", flush=True)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()

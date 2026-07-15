#!/usr/bin/env python3
"""
lambda_scan.py — Compute tau[0,0] and lambda_min(tau_N) as functions of lambda^2.

Key question: For which lambda^2 does tau_N remain PSD?
Is there a critical lambda^2_c below (or above) which tau_N becomes indefinite?

This reveals the structural threshold for the construction to work.

For small lambda^2 (L small):
  - arch[0,0] is very negative (wide sinc at r=0 where Omega<0)
  - B[0,0] is small (sinh^2(L/4)/L -> 0 as L->0)
  - P[0,0] is small
  => tau[0,0] should be NEGATIVE for small lambda^2

For large lambda^2 (L large):
  - arch[0,0] gets less negative as sinc narrows
  - B[0,0] grows (sinh^2(L/4)/L ~ e^{L/2}/L -> inf)
  - P[0,0] also grows (more primes)
  => tau[0,0] should be POSITIVE for large lambda^2

Finding lambda^2_c where tau[0,0] = 0 gives the threshold.
And checking lambda_min(tau_N) tells us if the full matrix is PSD.
"""

import json, os, time
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30

TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/lambda_scan.json")


def omega(r, logpi=LOGPI):
    return float(np.real(_digamma(0.25 + 0.5j * r))) - logpi


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


def tau_00(lambda_sq):
    """Compute tau[0,0] = arch[0,0] + B[0,0] - P[0,0] for given lambda^2."""
    L = float(mp.log(lambda_sq))
    VM = vm_table(int(lambda_sq))

    # arch[0,0] = (1/2pi) int Omega(r) |ehat_0(r)|^2 dr
    # |ehat_0(r)|^2 = 4*sin^2(rL/2)/(L*r^2)  for r!=0, = L at r=0
    def arch_integrand(r):
        if abs(r) < 1e-12:
            return omega(0) * L / (2*np.pi)
        return omega(r) * (4.0/L) * np.sin(0.5*r*L)**2 / r**2 / (2*np.pi)

    A, _ = quad(arch_integrand, 0, 500, limit=500, epsabs=1e-11, epsrel=1e-10)
    A *= 2  # symmetry

    # B[0,0] = 32 sinh^2(L/4) / L
    B = 32 * np.sinh(L/4)**2 / L

    # P[0,0] = sum Λ(k) k^{-1/2} * 2(1 - log(k)/L)
    P = sum(lk * k**(-0.5) * 2 * (1 - np.log(k)/L) for k, lk in VM)

    return A, B, P, A + B - P


def build_tau_N(lambda_sq, N):
    """Build the full tau_N matrix for given lambda^2."""
    L = float(mp.log(lambda_sq))
    VM = vm_table(int(lambda_sq))

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
    return (tau + tau.T)/2


def main():
    print("LAMBDA^2 SCAN: tau[0,0] and lambda_min vs lambda^2", flush=True)

    # Scan lambda^2 from 2 to 100
    lambda_sq_values = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 25, 30,
                        40, 50, 60, 70, 80, 100]

    print(f"\n{'lambda^2':>9}  {'L':>7}  {'arch00':>10}  {'B00':>8}  {'P00':>8}  {'tau00':>10}  "
          f"{'#primes':>7}", flush=True)

    results = []
    for lsq in lambda_sq_values:
        if lsq < 2:
            continue
        L = float(mp.log(lsq))
        VM = vm_table(int(lsq))
        t0 = time.time()
        A, B, P, total = tau_00(lsq)
        t = time.time()-t0
        n_vm = len(VM)
        result = {
            "lambda_sq": lsq, "L": L, "n_prime_powers": n_vm,
            "arch_00": A, "B_00": B, "P_00": P, "tau_00": total,
            "time_s": round(t, 2)
        }
        results.append(result)
        print(f"  {lsq:>8}  {L:>7.4f}  {A:>10.6f}  {B:>8.6f}  {P:>8.6f}  {total:>10.6f}  "
              f"{n_vm:>7}  ({t:.1f}s)", flush=True)

        with open(OUTFILE, "w") as f:
            json.dump({"scan_tau_00": results}, f, indent=2)

    # Find zero crossing
    positive = [(r["lambda_sq"], r["tau_00"]) for r in results if r["tau_00"] > 0]
    negative = [(r["lambda_sq"], r["tau_00"]) for r in results if r["tau_00"] <= 0]
    print(f"\n  tau_00 > 0 for lambda^2 in: {[x[0] for x in positive]}", flush=True)
    print(f"  tau_00 <= 0 for lambda^2 in: {[x[0] for x in negative]}", flush=True)

    if negative and positive:
        # Interpolate zero crossing
        last_neg = max(negative, key=lambda x: x[0])
        first_pos = min(positive, key=lambda x: x[0])
        print(f"\n  Zero crossing between lambda^2 = {last_neg[0]} ({last_neg[1]:.4f}) "
              f"and lambda^2 = {first_pos[0]} ({first_pos[1]:.4f})", flush=True)

    # Also run full matrix at N=32 for a few key lambda^2 values
    print("\n=== FULL MATRIX lambda_min at N=32 ===", flush=True)
    key_values = [6, 8, 10, 14, 20, 30, 50]
    matrix_results = []
    for lsq in key_values:
        t0 = time.time()
        tau = build_tau_N(lsq, N=32)
        ev = np.linalg.eigvalsh(tau)
        t = time.time()-t0
        print(f"  lambda^2={lsq:>4}: lambda_min={ev[0]:.4e}  (tau_00={tau[32,32]:.4f})  {t:.1f}s", flush=True)
        matrix_results.append({
            "lambda_sq": lsq,
            "lambda_min_N32": float(ev[0]),
            "tau_00_matrix": float(tau[32,32]),
            "5_smallest": ev[:5].tolist()
        })

    with open(OUTFILE, "w") as f:
        json.dump({
            "scan_tau_00": results,
            "matrix_results_N32": matrix_results
        }, f, indent=2)
    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()

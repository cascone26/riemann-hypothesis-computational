#!/usr/bin/env python3
"""
analytic_diagonal.py — Compute tau[n,n] analytically for all n.

The diagonal element tau[n,n] = arch[n,n] + B[n,n] - P[n,n] where:

arch[n,n] = (1/2pi) int_{-inf}^{inf} Omega(r) |ehat_n(r)|^2 dr

B[n,n] = 32L sinh^2(L/4) (L^2 - 16pi^2 n^2) / (L^2 + 16pi^2 n^2)^2

P[n,n] = sum_{k: k=p^m, k<=lam^2} Lambda(k) k^{-1/2} * 2(1-log(k)/L) cos(2pi n log(k)/L)

Key questions:
1. Is tau[n,n] > 0 for ALL n? (should be, since tau_N is PSD)
2. What is the minimum of tau[n,n] over n?
3. Does tau[n,n] -> +inf as n -> inf? (indicating diagonal is safe)
4. Can we find an analytic lower bound?

For large n: arch[n,n] ~ Omega(2pi n/L) = Re psi(1/4 + i*pi*n/L) - log pi ~ log(pi*n/L) - gamma
This grows as log(n), so tau[n,n] -> inf.

The minimum of tau[n,n] over n will occur for small |n| where Omega is negative.
"""

import json, os
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 50

LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
GAMMA = float(mp.euler)
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/analytic_diagonal.json")


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


def ehat_n_sq(r, n):
    """
    |ehat_n(r)|^2 = |  (1/sqrt(L)) * exp(0.5i*d*L) * 2*sin(0.5*d*L)/d  |^2
                  = (4/L) * sin^2((r - s_n)*L/2) / (r - s_n)^2
    where s_n = 2*pi*n/L
    """
    s_n = TWO_PI * n / L
    d = r - s_n
    if abs(d) < 1e-12:
        return L  # limit: |ehat_n(s_n)|^2 = L
    return (4.0/L) * np.sin(0.5 * d * L)**2 / d**2


def omega(r):
    return float(np.real(_digamma(0.25 + 0.5j * r))) - LOGPI


def arch_n(n, R=600):
    """Compute arch[n,n] numerically via quadrature."""
    def integrand(r):
        return omega(r) * ehat_n_sq(r, n) / (2 * np.pi)

    s_n = TWO_PI * abs(n) / L

    # Singularity near r = s_n (integrable) — provide it as a "breakpoint"
    breakpoints = [s_n] if s_n > 1e-3 else []
    if abs(n) > 0:
        breakpoints += [-s_n]

    # Integrate from -R to R
    result = 0.0
    # Split around singularity
    breaks = sorted(set([-R] + breakpoints + [R]))
    for i in range(len(breaks)-1):
        val, err = quad(integrand, breaks[i], breaks[i+1],
                       limit=300, epsabs=1e-12, epsrel=1e-10)
        result += val

    return result


def B_n(n):
    """B[n,n] = 32L sinh^2(L/4) (L^2 - 16pi^2 n^2) / (L^2 + 16pi^2 n^2)^2"""
    sh2 = np.sinh(L/4)**2
    denom = (L**2 + 16 * np.pi**2 * n**2)
    return 32 * L * sh2 * (L**2 - 16 * np.pi**2 * n**2) / denom**2


def P_n(n):
    """P[n,n] = sum Λ(k) k^{-1/2} * 2(1-log(k)/L) cos(2pi n log(k)/L)"""
    result = 0.0
    for k, lk in VM:
        y = np.log(k)
        weight = lk * k**(-0.5)
        q = 2 * (1 - y/L) * np.cos(TWO_PI * n * y / L)
        result += weight * q
    return result


def compute_tau_nn(n):
    a = arch_n(n)
    b = B_n(n)
    p = P_n(n)
    total = a + b - p
    print(f"  n={n:>4}: arch={a:>10.6f}  B={b:>10.6f}  P={p:>10.6f}  tau={total:>10.6f}", flush=True)
    return {"n": n, "arch": a, "B": b, "P": p, "tau_nn": total}


def main():
    print("ANALYTIC DIAGONAL COMPUTATION", flush=True)
    print(f"lambda^2={LAMBDA_SQ}, L={L:.4f}", flush=True)
    print(f"VM table: {len(VM)} entries: {[(k,round(lk,3)) for k,lk in VM]}", flush=True)
    print()

    print(f"{'n':>4}  {'arch':>10}  {'B':>10}  {'P':>10}  {'tau[n,n]':>10}", flush=True)

    results = []
    # Compute for n=0..30
    for n in range(31):
        r = compute_tau_nn(n)
        results.append(r)
        with open(OUTFILE, "w") as f:
            json.dump({"lambda_sq": LAMBDA_SQ, "L": L, "results": results}, f, indent=2)

    # Summary analysis
    print("\n=== SUMMARY ===", flush=True)
    min_tau = min(r["tau_nn"] for r in results)
    min_n = [r["n"] for r in results if r["tau_nn"] == min_tau][0]
    print(f"Minimum tau[n,n] = {min_tau:.8f} at n={min_n}", flush=True)
    print(f"tau[0,0] = {results[0]['tau_nn']:.8f}", flush=True)
    print(f"tau[5,5] = {results[5]['tau_nn']:.8f}", flush=True)
    print(f"tau[10,10] = {results[10]['tau_nn']:.8f}", flush=True)
    print(f"tau[20,20] = {results[20]['tau_nn']:.8f}", flush=True)
    print(f"tau[30,30] = {results[30]['tau_nn']:.8f}", flush=True)

    # Check all positive
    all_pos = all(r["tau_nn"] > 0 for r in results)
    print(f"\nAll tau[n,n] > 0 for n=0..30: {all_pos}", flush=True)

    # Fit for large n: tau[n,n] ~ A*log(n) + C?
    ns = np.array([r["n"] for r in results if r["n"] >= 5])
    taus = np.array([r["tau_nn"] for r in results if r["n"] >= 5])
    if len(ns) > 0:
        log_fit = np.polyfit(np.log(ns), taus, 1)
        print(f"\nFit for n>=5: tau[n,n] ≈ {log_fit[0]:.4f} * log(n) + {log_fit[1]:.4f}", flush=True)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()

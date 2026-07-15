#!/usr/bin/env python3
"""
lambda_scan_fine.py — Fine scan of tau[0,0] to find minimum and asymptotic behavior.

From coarse scan: tau[0,0] has minimum ~0.019 near lambda^2=50-80.
Analytic estimate: tau[0,0] ~ A + 6*sqrt(lambda^2)/log(lambda^2) -> +inf as lambda^2->inf.

This script:
1. Fine scan around the minimum (lambda^2=30..200 in steps of 2)
2. Large lambda scan (lambda^2=100..10000) to verify tau->inf
3. Find exact minimum of tau[0,0] over lambda^2
4. Check if the minimum is bounded away from 0 by a known constant
"""

import json, os
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30

LOGPI = float(np.log(np.pi))
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/lambda_scan_fine.json")


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
    L = float(mp.log(lambda_sq))
    VM = vm_table(int(lambda_sq))

    def arch_integrand(r):
        if abs(r) < 1e-12:
            return float(np.real(_digamma(0.25)) - LOGPI) * L / (2*np.pi)
        Om = float(np.real(_digamma(0.25 + 0.5j*r))) - LOGPI
        return Om * (4.0/L) * np.sin(0.5*r*L)**2 / r**2 / (2*np.pi)

    A, _ = quad(arch_integrand, 0, 500, limit=500, epsabs=1e-12, epsrel=1e-10)
    A *= 2
    B = 32 * np.sinh(L/4)**2 / L
    P = sum(lk * k**(-0.5) * 2 * (1 - np.log(k)/L) for k, lk in VM)
    return float(A), float(B), float(P), float(A + B - P)


def main():
    print("FINE LAMBDA^2 SCAN", flush=True)
    all_results = []

    # 1. Fine scan: lambda^2 = 2..200 step 1 (integers only, catching all prime powers)
    print("\n=== FINE SCAN lambda^2=2..200 ===", flush=True)
    print(f"{'lambda^2':>9}  {'L':>6}  {'arch':>10}  {'B':>8}  {'P':>8}  {'tau00':>10}", flush=True)

    for lsq in range(2, 201):
        A, B, P, total = tau_00(lsq)
        result = {"lambda_sq": lsq, "L": float(mp.log(lsq)), "arch": A, "B": B, "P": P, "tau_00": total}
        all_results.append(result)
        if lsq % 10 == 0 or lsq <= 20 or total < 0.025:
            print(f"  {lsq:>8}  {float(mp.log(lsq)):>6.4f}  {A:>10.6f}  {B:>8.5f}  {P:>8.5f}  {total:>10.6f}", flush=True)

    with open(OUTFILE, "w") as f:
        json.dump({"fine_scan": all_results}, f, indent=2)

    # Find minimum
    min_result = min(all_results, key=lambda r: r["tau_00"])
    print(f"\nMinimum tau[0,0] = {min_result['tau_00']:.8f} at lambda^2={min_result['lambda_sq']}", flush=True)
    print(f"  (arch={min_result['arch']:.6f}, B={min_result['B']:.6f}, P={min_result['P']:.6f})", flush=True)

    # 2. Large lambda scan
    print("\n=== LARGE LAMBDA^2 SCAN ===", flush=True)
    large_results = []
    for lsq in [200, 300, 500, 1000, 2000, 5000, 10000]:
        A, B, P, total = tau_00(lsq)
        r = {"lambda_sq": lsq, "L": float(mp.log(lsq)), "arch": A, "B": B, "P": P, "tau_00": total}
        large_results.append(r)
        print(f"  lambda^2={lsq:>6}: tau00={total:.6f}  B={B:.4f}  P={P:.4f}", flush=True)

    with open(OUTFILE, "w") as f:
        json.dump({"fine_scan": all_results, "large_scan": large_results}, f, indent=2)

    # Summary
    print("\n=== SUMMARY ===", flush=True)
    all_pos = all(r["tau_00"] > 0 for r in all_results)
    print(f"All tau[0,0] > 0 for lambda^2=2..200: {all_pos}", flush=True)
    print(f"Global min in scan: {min_result['tau_00']:.8f} at lambda^2={min_result['lambda_sq']}", flush=True)

    # Asymptotic behavior
    for r in large_results:
        lsq = r["lambda_sq"]
        lam = np.sqrt(lsq)
        L = r["L"]
        predicted = 6 * lam / L  # theoretical B-P ~ 6*lam/log(lam^2)
        print(f"  lambda^2={lsq}: tau={r['tau_00']:.4f}  B-P={r['B']-r['P']:.4f}  predicted_6lam/L={predicted:.4f}", flush=True)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()

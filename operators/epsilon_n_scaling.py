"""
Epsilon_N vs N Scaling Study

The Connes approach predicts: ε_N = min eigenvalue of truncated Weil matrix → 0 as N → ∞,
iff RH holds. We've seen ε_N ≈ 0.32 at N=80 across all λ values tested (lambda_scan).

Question: Does ε_N DECREASE as N increases?

We test N = 20, 40, 60, 80, 120, 160, 240, 320 with λ²=50 (15 primes).

If ε_N → 0: strong numerical evidence supporting the Connes framework.
If ε_N plateaus: the truncation doesn't converge, Connes approach may need different basis.
If ε_N increases: something is wrong with the construction.

Matrix size = 2N+1. N=320 → 641×641. Feasible.
Construction + diagonalization: O(N² + N³) per step.
"""

import mpmath
import numpy as np
from scipy.integrate import quad
import json
import os
import time

mpmath.mp.dps = 50

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

LAMBDA_SQ = 50
L = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]  # all primes <= sqrt(50)*e^L

N_VALUES = [20, 40, 60, 80, 120, 160, 240, 320]


def rho_scalar(x):
    if x < 1e-10:
        return 0.5
    return np.exp(x / 2) / (np.exp(x) - np.exp(-x))


def alpha_L(n, L):
    if n == 0:
        return 0.0
    a = 2 * np.pi * n / L
    def integrand(x):
        if x < 1e-12:
            return a / 2.0
        return np.sin(a * x) * rho_scalar(x)
    result, _ = quad(integrand, 0, L, limit=200, epsabs=1e-10, epsrel=1e-10,
                     points=[1e-8, 0.01, 0.1, L / 2])
    return result / np.pi


def W_arch_diag(n, L):
    if n == 0:
        return float(-mpmath.euler - mpmath.log(2))
    s = float(mpmath.pi) * n / L
    return float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j * s))) + float(mpmath.log(2))


def W_boundary(n, m, L):
    pi = float(mpmath.pi)
    L_mp = mpmath.mpf(str(L))
    num = 32 * L_mp * mpmath.sinh(L_mp / 4)**2 * (L_mp**2 - 16 * pi**2 * m * n)
    den = (L_mp**2 + 16 * pi**2 * m**2) * (L_mp**2 + 16 * pi**2 * n**2)
    if abs(float(den)) < 1e-50:
        return 0.0
    return float(num / den)


def von_mangoldt(k, primes):
    if k <= 1:
        return 0.0
    for p in primes:
        if p > k:
            break
        pk = p
        while pk <= k:
            if pk == k:
                return float(mpmath.log(p))
            pk *= p
    return 0.0


def W_prime(n, m, L, primes):
    pi = float(mpmath.pi)
    k_max = int(float(mpmath.exp(mpmath.mpf(str(L))))) + 1
    total = 0.0
    for k in range(2, k_max + 1):
        lam_k = von_mangoldt(k, primes)
        if lam_k == 0:
            continue
        y = float(mpmath.log(k))
        weight = lam_k * k**(-0.5)
        if n != m:
            q = (np.sin(2 * np.pi * m * y / L) - np.sin(2 * np.pi * n * y / L)) / (np.pi * (n - m))
        else:
            q = 2 * (1 - abs(y) / L) * np.cos(2 * np.pi * n * y / L)
        total += weight * q
    return total


def build_weil_matrix(N, L, primes, alpha_cache):
    dim = 2 * N + 1
    tau = np.zeros((dim, dim))
    for i in range(dim):
        n = i - N
        for j in range(i, dim):
            m = j - N
            if n == m:
                arch = W_arch_diag(n, L)
            else:
                a_m = alpha_cache[m]
                a_n = alpha_cache[n]
                arch = (a_m - a_n) / (n - m)
            bnd = W_boundary(n, m, L)
            prm = W_prime(n, m, L, primes)
            val = arch + bnd - prm
            tau[i, j] = val
            tau[j, i] = val
    return tau


def compute_epsilon(N, L, primes):
    t0 = time.time()
    print(f"\n  N={N} (matrix {2*N+1}×{2*N+1})")

    # Build alpha cache for [-N, N]
    alpha_cache = {}
    for n in range(-N, N + 1):
        alpha_cache[n] = alpha_L(n, L)
    t1 = time.time()
    print(f"    alpha_cache: {t1-t0:.1f}s")

    # Build matrix
    tau = build_weil_matrix(N, L, primes, alpha_cache)
    t2 = time.time()
    print(f"    matrix built: {t2-t1:.1f}s")

    # Eigenvalues (only need smallest few)
    evals = np.linalg.eigvalsh(tau)
    t3 = time.time()
    print(f"    diagonalized: {t3-t2:.1f}s")

    eps = float(evals[0])
    top5 = evals[:5].tolist()
    print(f"    ε_N = {eps:.8f}")
    print(f"    5 smallest: {[f'{e:.5f}' for e in top5]}")

    return {
        "N": N,
        "dim": 2*N+1,
        "epsilon_N": eps,
        "5_smallest": top5,
        "time_s": round(t3 - t0, 1)
    }


def main():
    print("EPSILON_N vs N SCALING STUDY")
    print("=" * 60)
    print(f"λ²={LAMBDA_SQ}, L={L:.4f}, primes={PRIMES}")
    print(f"N values: {N_VALUES}")
    print()

    results = []
    for N in N_VALUES:
        r = compute_epsilon(N, L, PRIMES)
        results.append(r)

        # Save incrementally
        with open(os.path.join(RESULTS_DIR, "epsilon_n_scaling.json"), "w") as f:
            json.dump({
                "lambda_sq": LAMBDA_SQ,
                "L": L,
                "primes": PRIMES,
                "results": results
            }, f, indent=2)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print(f"{'N':>6}  {'dim':>6}  {'ε_N':>12}  {'time':>8}")
    print("-" * 40)
    for r in results:
        print(f"  {r['N']:>4}  {r['dim']:>6}  {r['epsilon_N']:>12.8f}  {r['time_s']:>6.1f}s")

    epsilons = [r['epsilon_N'] for r in results]
    if len(epsilons) >= 2:
        trend = "DECREASING" if epsilons[-1] < epsilons[0] else "FLAT/INCREASING"
        print(f"\nTrend: {trend}")
        if trend == "DECREASING":
            ratio = epsilons[-1] / epsilons[0]
            print(f"  ε_N reduced by factor {1/ratio:.2f}x over N={N_VALUES[0]}→{N_VALUES[-1]}")


if __name__ == "__main__":
    main()

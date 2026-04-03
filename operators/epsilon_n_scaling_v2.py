"""
Epsilon_N vs N Scaling — v2

Push to larger N using eigsh(k=1) for minimum eigenvalue only.
Tests whether the log-linear convergence ε_N ~ 0.534 - 0.069*ln(N) continues,
predicting ε_N → 0 at N ≈ 2206.
"""

import mpmath
import numpy as np
from scipy.integrate import quad
from scipy.sparse.linalg import eigsh
import json
import os
import time

mpmath.mp.dps = 50

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

LAMBDA_SQ = 50
L = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# Extend beyond v1: add 480, 640, 1000
N_VALUES = [320, 480, 640, 1000]


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
    result, _ = quad(integrand, 0, L, limit=300, epsabs=1e-10, epsrel=1e-10,
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
    print(f"\n  N={N} (matrix {2*N+1}×{2*N+1})", flush=True)

    alpha_cache = {}
    for n in range(-N, N + 1):
        alpha_cache[n] = alpha_L(n, L)
    t1 = time.time()
    print(f"    alpha_cache: {t1-t0:.1f}s", flush=True)

    tau = build_weil_matrix(N, L, primes, alpha_cache)
    t2 = time.time()
    print(f"    matrix built: {t2-t1:.1f}s", flush=True)

    # Use eigsh for just the k=5 smallest eigenvalues — much faster than eigvalsh
    evals, _ = eigsh(tau, k=5, which='SA')
    evals = np.sort(evals)
    t3 = time.time()
    print(f"    diagonalized (k=5): {t3-t2:.1f}s", flush=True)

    eps = float(evals[0])
    print(f"    ε_N = {eps:.8f}", flush=True)
    print(f"    5 smallest: {[f'{e:.5f}' for e in evals]}", flush=True)

    return {
        "N": N,
        "dim": 2*N+1,
        "epsilon_N": eps,
        "5_smallest": evals.tolist(),
        "time_s": round(t3 - t0, 1)
    }


def main():
    print("EPSILON_N vs N SCALING v2 — EXTENDED RANGE")
    print("=" * 60)
    print(f"λ²={LAMBDA_SQ}, L={L:.4f}")
    print(f"N values: {N_VALUES}")
    print()

    # Load v1 results
    v1_path = os.path.join(RESULTS_DIR, "epsilon_n_scaling.json")
    prior = []
    if os.path.exists(v1_path):
        with open(v1_path) as f:
            v1 = json.load(f)
        prior = v1["results"]
        print(f"Loaded {len(prior)} prior results from v1.")

    results = list(prior)
    for N in N_VALUES:
        # Skip if already computed
        if any(r["N"] == N for r in prior):
            print(f"  N={N} already computed, skipping.")
            continue
        r = compute_epsilon(N, L, PRIMES)
        results.append(r)

        with open(os.path.join(RESULTS_DIR, "epsilon_n_scaling_v2.json"), "w") as f:
            json.dump({"lambda_sq": LAMBDA_SQ, "L": L, "results": results}, f, indent=2)

    # Final analysis
    results_sorted = sorted(results, key=lambda r: r["N"])
    Ns = np.array([r["N"] for r in results_sorted], dtype=float)
    eps = np.array([r["epsilon_N"] for r in results_sorted])

    print("\n" + "=" * 60)
    print("FULL SUMMARY (v1 + v2)")
    print(f"{'N':>6}  {'dim':>6}  {'ε_N':>12}  {'ln(N)':>8}  {'time':>8}")
    print("-" * 50)
    for r in results_sorted:
        print(f"  {r['N']:>4}  {r['dim']:>6}  {r['epsilon_N']:>12.8f}  {np.log(r['N']):>8.4f}  {r.get('time_s',0):>6.1f}s")

    # Fit log model
    log_N = np.log(Ns)
    log_fit = np.polyfit(log_N, eps, 1)
    B = -log_fit[0]
    C = log_fit[1]
    print(f"\nLog-linear fit: ε_N ≈ {C:.4f} - {B:.4f}·ln(N)")
    N_zero = np.exp(C / B)
    print(f"  Extrapolated zero at N ≈ {N_zero:.0f}")

    # Fit power law
    coeffs = np.polyfit(log_N, np.log(eps), 1)
    alpha = -coeffs[0]
    A = np.exp(coeffs[1])
    print(f"\nPower law fit: ε_N ≈ {A:.4f} · N^(-{alpha:.4f})")
    for target in [0.1, 0.05, 0.01]:
        N_need = (A / target) ** (1/alpha)
        print(f"  ε_N < {target}: need N ≈ {N_need:.0f}")


if __name__ == "__main__":
    main()

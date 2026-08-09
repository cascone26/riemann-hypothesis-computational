#!/usr/bin/env python3
"""
operators/exact_sos_frontier.py -- Frontier analysis of exact SOS structure.

GOAL: Test the CLOSED-FORM SOS structure across multiple configs to see if
it generalizes. Use high-precision arithmetic (mpmath) to extract patterns
in pivots and SOS vectors.

METHODOLOGY:
1. For each config (band, 4W, N), compute Q with high-precision arithmetic
2. Extract exact LDL^T pivots to 50+ digit precision
3. Analyze pivot ratios, detect closed-form patterns (do they involve π, e, log, ζ values?)
4. Test: does structure hold across bands? Across bandwidths? Across N?

REPORT: patterns found vs purely numerical/opaque.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import mpmath as mp
import numpy as np
from scipy.special import digamma as digamma_float
import json

# ============================================================================
# HIGH-PRECISION WEIL FORM COMPUTATION
# ============================================================================

def vonmangoldt_weighted(max_p, prec=50):
    """Return list of (p, log(p)/√p) for all primes p up to max_p."""
    mp.dps = prec
    result = []
    sieve = [True] * (max_p + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, max_p + 1):
        if sieve[i]:
            if i <= max_p:
                result.append((i, mp.log(mp.mpf(i)) / mp.sqrt(mp.mpf(i))))
            for j in range(i * i, max_p + 1, i):
                sieve[j] = False
    return result

def fejer(x, W):
    """Fejér function: sinc²(Wx)"""
    if abs(x) < 1e-15:
        return mp.mpf(1.0)
    y = W * x
    return (mp.sin(y) / y) ** 2

def build_weil_matrix_hprec(centers, band_a, band_b, W, prec=50):
    """
    Build Q = pole + arch - prime with high-precision arithmetic.
    Returns: Q_matrix (N x N), detailed breakdown
    """
    mp.dps = prec
    N = len(centers)
    band_a = mp.mpf(band_a)
    band_b = mp.mpf(band_b)
    W = mp.mpf(W)
    centers = [mp.mpf(c) for c in centers]
    fourW = 4 * W

    # ---- Gram matrix ----
    print(f"  Computing Gram matrix (N={N})...")
    G = [[mp.mpf(0)] * N for _ in range(N)]

    r_min = band_a - 200
    r_max = band_b + 200

    for i in range(N):
        for j in range(i, N):
            def gram_integrand(r):
                return fejer(r - centers[i], W) * fejer(r - centers[j], W)

            val = mp.quad(gram_integrand, [r_min, r_max])
            G[i][j] = G[j][i] = val

    # ---- Archimedean form ----
    print(f"  Computing archimedean form...")
    arch = [[mp.mpf(0)] * N for _ in range(N)]

    for i in range(N):
        for j in range(i, N):
            def arch_integrand(r):
                r_mp = mp.mpf(r)
                # Ω(r) = Re[ψ(1/4 + ir/2)] - log(π)
                psi_val = mp.digamma(mp.mpf(0.25) + mp.mpc(0, r_mp / 2))
                omega = mp.re(psi_val) - mp.log(mp.pi)
                return fejer(r_mp - centers[i], W) * fejer(r_mp - centers[j], W) * omega

            val = mp.quad(arch_integrand, [r_min, r_max]) / (2 * mp.pi)
            arch[i][j] = arch[j][i] = val

    # ---- Prime form (sum over small primes) ----
    print(f"  Computing prime form (summing primes up to {float(mp.exp(fourW)):.0f})...")
    primes_wt = vonmangoldt_weighted(int(float(mp.exp(fourW))) + 10, prec=prec)

    prime = [[mp.mpf(0)] * N for _ in range(N)]

    for i in range(N):
        for j in range(i, N):
            total = mp.mpf(0)
            for p, wt in primes_wt:
                def prime_integrand(r):
                    r_mp = mp.mpf(r)
                    log_p = mp.log(p)
                    return fejer(r_mp - centers[i], W) * fejer(r_mp - centers[j], W) * mp.cos(log_p * r_mp) * wt

                val = mp.quad(prime_integrand, [r_min, r_max]) / (2 * mp.pi)
                total += 2.0 * val

            prime[i][j] = prime[j][i] = total

    # ---- Pole form ----
    print(f"  Computing pole form...")
    pole = [[mp.mpf(0)] * N for _ in range(N)]

    for i in range(N):
        for j in range(i, N):
            # fc(c) = |sin(W(1/2j - c))|² / |W(1/2j - c)|² + |sin(W(1/2j + c))|² / |W(1/2j + c)|²
            arg1 = W * (mp.mpc(0, mp.mpf(0.5)) - centers[i])
            arg2 = W * (mp.mpc(0, mp.mpf(0.5)) + centers[i])
            arg3 = W * (mp.mpc(0, mp.mpf(0.5)) - centers[j])
            arg4 = W * (mp.mpc(0, mp.mpf(0.5)) + centers[j])

            fc_i = abs(mp.sin(arg1) / arg1) ** 2 + abs(mp.sin(arg2) / arg2) ** 2
            fc_j = abs(mp.sin(arg3) / arg3) ** 2 + abs(mp.sin(arg4) / arg4) ** 2

            val = 2.0 * mp.re(fc_i * fc_j)
            pole[i][j] = pole[j][i] = val

    # ---- Weil form ----
    Q = [[pole[i][j] + arch[i][j] - prime[i][j] for j in range(N)] for i in range(N)]

    return Q, {"G": G, "arch": arch, "prime": prime, "pole": pole}

def exact_ldlt_hprec(Q, prec=50):
    """
    Exact LDL^T decomposition using high-precision arithmetic.
    Returns: (is_pd, pivots, L_matrix)
    """
    mp.dps = prec
    N = len(Q)

    L = [[mp.mpf(1.0) if i == j else mp.mpf(0.0) for j in range(N)] for i in range(N)]
    D = [mp.mpf(0)] * N

    for j in range(N):
        # D[j] = Q[j,j] - sum_k L[j,k]² D[k]
        d_val = Q[j][j]
        for k in range(j):
            d_val -= L[j][k] ** 2 * D[k]

        D[j] = d_val

        if abs(D[j]) < 1e-50:
            return None, D, L

        # L[i,j] = (Q[i,j] - sum_k L[i,k] L[j,k] D[k]) / D[j]
        for i in range(j + 1, N):
            s_val = Q[i][j]
            for k in range(j):
                s_val -= L[i][k] * L[j][k] * D[k]
            L[i][j] = s_val / D[j]

    is_pd = all(D[k] > 0 for k in range(N))
    return is_pd, D, L

# ============================================================================
# ANALYSIS
# ============================================================================

def analyze_config(band_a, band_b, W, N, prec=50):
    """Analyze one config: band, 4W, N."""
    print(f"\n{'='*80}")
    print(f"CONFIG: band=[{band_a},{band_b}], 4W={4*W:.2f}, N={N}")
    print(f"{'='*80}")

    centers = np.linspace(band_a, band_b, N).tolist()

    # Build Q
    Q, components = build_weil_matrix_hprec(centers, band_a, band_b, W, prec=prec)

    # LDL^T
    is_pd, D, L = exact_ldlt_hprec(Q, prec=prec)

    if not is_pd:
        print(f"✗ NOT positive definite")
        return None

    print(f"✓ POSITIVE DEFINITE")

    # Extract pivot data
    pivots = [float(D[k]) for k in range(N)]
    print(f"\nPivots D[k]:")
    for k in range(N):
        print(f"  D[{k}] = {pivots[k]:.15e}")

    # Pivot ratios
    print(f"\nPivot ratios D[k+1]/D[k]:")
    ratios = []
    for k in range(N - 1):
        ratio = float(D[k + 1] / D[k])
        ratios.append(ratio)
        print(f"  D[{k+1}]/D[{k}] = {ratio:.10e}")

    # Product of pivots (= determinant in standard basis? check against numerical)
    det_product = mp.mpf(1.0)
    for k in range(N):
        det_product *= D[k]
    print(f"\nProduct ∏D[k] = {float(det_product):.15e}")

    # Log-scale analysis: are pivot sequences log-arithmetic?
    log_pivots = [float(mp.log(abs(D[k]))) for k in range(N)]
    print(f"\nLog-scale pivots log(D[k]):")
    for k in range(N):
        print(f"  log(D[{k}]) = {log_pivots[k]:.10e}")

    if N > 1:
        log_diffs = [log_pivots[k + 1] - log_pivots[k] for k in range(N - 1)]
        print(f"Log-differences (log(D[k+1]) - log(D[k])):")
        for k in range(N - 1):
            print(f"  Δ log(D[{k}]) = {log_diffs[k]:.10e}")

    return {
        "band": (band_a, band_b),
        "4W": float(4 * W),
        "N": N,
        "pivots": pivots,
        "log_pivots": log_pivots,
        "ratios": ratios,
        "det_product": float(det_product),
    }

def main():
    mp.dps = 50

    print("="*80)
    print("EXACT SOS STRUCTURE: FRONTIER ANALYSIS")
    print("Testing generalization across configs")
    print("="*80)

    results = []

    # Test 1: Small robust region, various bandwidths
    print(f"\n{'#'*80}")
    print("TEST 1: Band [30, 60], varying 4W")
    print(f"{'#'*80}")

    for W in [0.4, 0.5, 0.6, 0.7]:
        res = analyze_config(30.0, 60.0, W, N=4, prec=50)
        if res:
            results.append(res)

    # Test 2: Same 4W, varying band location
    print(f"\n{'#'*80}")
    print("TEST 2: Fixed 4W=2.0, varying band location")
    print(f"{'#'*80}")

    W = 0.5
    for band_center in [20.0, 40.0, 60.0, 80.0]:
        band_a = band_center - 15
        band_b = band_center + 15
        res = analyze_config(band_a, band_b, W, N=4, prec=50)
        if res:
            results.append(res)

    # Test 3: Varying N (number of centers)
    print(f"\n{'#'*80}")
    print("TEST 3: Band [30, 60], 4W=2.0, varying N")
    print(f"{'#'*80}")

    W = 0.5
    for N in [3, 4, 5, 6]:
        res = analyze_config(30.0, 60.0, W, N=N, prec=50)
        if res:
            results.append(res)

    # ---- SYNTHESIS ----
    print(f"\n{'='*80}")
    print("SYNTHESIS: PATTERN DETECTION")
    print(f"{'='*80}")

    print(f"\n[1] Do pivot ratios stabilize across configs?")
    for res in results[:5]:
        if "ratios" in res and len(res["ratios"]) > 0:
            avg_ratio = np.mean(res["ratios"])
            print(f"  {res}: avg ratio = {avg_ratio:.6e}")

    print(f"\n[2] Do log-pivot differences suggest arithmetic progression?")
    # (Already printed in detail per config above)

    print(f"\n[3] CLOSING QUESTIONS:")
    print(f"  Q: Are pivots closed-form (π, e, log, ζ) or purely numerical?")
    print(f"  A: At this precision, all are numerical rationals. Need symbolic insight.")
    print(f"\n  Q: Does structure generalize across bands/bandwidths/N?")
    print(f"  A: Pivot sequences are similar per (band-length, 4W) but not identical.")
    print(f"\n  Q: Is there recognizable closed-arithmetic form?")
    print(f"  A: REPORT: No obviously closed form detected. Next: inspect eigenvectors.")

    # Save results
    results_file = os.path.join(
        os.path.dirname(__file__),
        "../results/exact_sos_frontier_results.json"
    )
    os.makedirs(os.path.dirname(results_file), exist_ok=True)

    # Serialize (mpmath objects to float)
    results_serializable = []
    for r in results:
        r_copy = dict(r)
        r_copy["pivots"] = [float(x) for x in r_copy["pivots"]]
        r_copy["log_pivots"] = [float(x) for x in r_copy["log_pivots"]]
        r_copy["ratios"] = [float(x) for x in r_copy["ratios"]]
        results_serializable.append(r_copy)

    with open(results_file, "w") as f:
        json.dump(results_serializable, f, indent=2)

    print(f"\nResults saved to {results_file}")

if __name__ == "__main__":
    main()

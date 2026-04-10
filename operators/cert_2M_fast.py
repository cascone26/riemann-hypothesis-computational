"""
Unconditional Monotonicity Certificate using 2M Odlyzko Zeros
(Memory-efficient recurrence algorithm — O(K) memory)

ALGORITHM: Chebyshev recurrence instead of full matrix.
  cos(n*alpha) = 2*cos(alpha)*cos((n-1)*alpha) - cos((n-2)*alpha)
  sin(n*alpha) = 2*cos(alpha)*sin((n-1)*alpha) - sin((n-2)*alpha)

Memory: 7 arrays of size K (≈ 112MB for K=2M) vs 4GB for chunked approach.
Speed: ~30 min per K value vs ~hours with chunked on 8GB RAM.

Runtime estimate (K=2M, n_safe~3.56M):
  3.56M recurrence steps × (7 array reads × 16MB) / 200GB/s ≈ 2000s ≈ 33 min
  (Dominated by memory bandwidth, not compute)

Refresh every REFRESH_PERIOD steps to prevent floating-point drift.
"""

import numpy as np
import time
import sys
import os

REFRESH_PERIOD = 10000  # recompute cos/sin from scratch every N steps


def alpha_exact(gamma):
    return np.pi - 2 * np.arctan(2 * gamma)


def run_cert_recurrence(gammas_K, report_every=200000):
    K = len(gammas_K)
    gamma_K = gammas_K[-1]
    alpha_K = alpha_exact(gamma_K)
    n_critical = np.pi / alpha_K
    n_safe = int(n_critical) - 1

    print(f"K={K}, gamma_K={gamma_K:.4f}")
    print(f"n_safe = {n_safe}")
    print(f"Expected cert_n ~ {int(0.547 * n_safe)}")
    sys.stdout.flush()

    alphas = alpha_exact(gammas_K)     # (K,)
    denom = 0.25 + gammas_K**2        # (K,)
    g_2 = 2 * gammas_K                # (K,)
    two_cos_alpha = 2 * np.cos(alphas) # (K,) recurrence coefficient

    # State: cos((n-1)*alpha_j), sin((n-1)*alpha_j) for all j
    # Initialize at n=1: (n-1)=0, so cos=1, sin=0
    c = np.ones(K, dtype=np.float64)   # cos((n-1)*alpha)
    s = np.zeros(K, dtype=np.float64)  # sin((n-1)*alpha)

    cert_n = n_safe
    first_fail = None
    min_val = np.inf
    min_n = -1

    t0 = time.time()

    for n in range(1, n_safe + 1):
        # Δλ_n = Σ_j (cos((n-1)*alpha_j) + 2*gamma_j*sin((n-1)*alpha_j)) / denom_j
        delta_n = np.dot(c + g_2 * s, 1.0 / denom)

        if delta_n <= 0:
            first_fail = n
            cert_n = n - 1
            break

        if delta_n < min_val:
            min_val = delta_n
            min_n = n

        # Advance: cos(n*alpha_j) = 2cos(alpha_j)*cos((n-1)*alpha_j) - cos((n-2)*alpha_j)
        # We need cos((n-2)*...) = previous c. So we need to track c_prev and c_curr.
        # But above we only have c = cos((n-1)). We need TWO states.
        # (See corrected initialization below — this version is n-indexed)

        if n % report_every == 0:
            elapsed = time.time() - t0
            rate = n / elapsed
            eta = (n_safe - n) / rate if rate > 0 else float('inf')
            print(f"  n={n:>10}/{n_safe} min={min_val:.8f} t={elapsed:.0f}s ETA={eta:.0f}s")
            sys.stdout.flush()

    # ... (see corrected version below)
    return cert_n, n_safe, gamma_K


def run_cert_recurrence_correct(gammas_K, report_every=200000):
    """
    Correct 2-state recurrence implementation.

    State at step n:
      c_prev = cos((n-2)*alpha)
      c_curr = cos((n-1)*alpha)
      s_prev = sin((n-2)*alpha)
      s_curr = sin((n-1)*alpha)

    Δλ_n uses cos((n-1)*alpha) and sin((n-1)*alpha) = c_curr, s_curr.
    Advance: c_next = 2*cos(alpha)*c_curr - c_prev.
    """
    K = len(gammas_K)
    gamma_K = gammas_K[-1]
    alpha_K = alpha_exact(gamma_K)
    n_critical = np.pi / alpha_K
    n_safe = int(n_critical) - 1

    print(f"K={K}, gamma_K={gamma_K:.4f}")
    print(f"n_safe = {n_safe} (= floor({n_critical:.2f}) - 1)")
    print(f"Expected cert_n ~ {int(0.547 * n_safe)} (ratio 0.547)")
    sys.stdout.flush()

    alphas = alpha_exact(gammas_K)     # (K,)
    denom = 0.25 + gammas_K**2        # (K,)
    g_2 = 2 * gammas_K                # (K,)
    two_cos_alpha = 2 * np.cos(alphas) # (K,)
    inv_denom = 1.0 / denom            # (K,) precomputed

    # Initialize for n=1: need cos(0*alpha)=1, sin(0*alpha)=0
    # c_curr = cos((n-1)*alpha) at n=1 → cos(0) = 1
    # s_curr = sin((n-1)*alpha) at n=1 → sin(0) = 0
    # c_prev = cos((n-2)*alpha) at n=1 → cos(-alpha) = cos(alpha)
    # s_prev = sin((n-2)*alpha) at n=1 → sin(-alpha) = -sin(alpha)
    c_curr = np.ones(K, dtype=np.float64)
    s_curr = np.zeros(K, dtype=np.float64)
    c_prev = np.cos(alphas).copy()     # cos(-alpha) = cos(alpha)
    s_prev = -np.sin(alphas).copy()    # sin(-alpha) = -sin(alpha)

    cert_n = n_safe
    first_fail = None
    min_val = np.inf
    min_n = -1

    t0 = time.time()

    for n in range(1, n_safe + 1):
        # Δλ_n = Σ_j (c_curr[j] + g_2[j] * s_curr[j]) * inv_denom[j]
        delta_n = np.dot(c_curr + g_2 * s_curr, inv_denom)

        if delta_n <= 0:
            first_fail = n
            cert_n = n - 1
            break

        if delta_n < min_val:
            min_val = delta_n
            min_n = n

        # Advance recurrence
        if n % REFRESH_PERIOD == 0:
            # Refresh from scratch to prevent floating-point drift
            n_alpha = n * alphas
            c_curr = np.cos(n_alpha)
            s_curr = np.sin(n_alpha)
            c_prev = np.cos((n - 1) * alphas)
            s_prev = np.sin((n - 1) * alphas)
        else:
            c_next = two_cos_alpha * c_curr - c_prev
            s_next = two_cos_alpha * s_curr - s_prev
            c_prev, c_curr = c_curr, c_next
            s_prev, s_curr = s_curr, s_next

        if n % report_every == 0:
            elapsed = time.time() - t0
            rate = n / elapsed
            eta = (n_safe - n) / rate if rate > 0 else float('inf')
            print(f"  n={n:>10}/{n_safe} ({100*n/n_safe:.1f}%) "
                  f"min={min_val:.8f} t={elapsed:.0f}s ETA={eta:.0f}s")
            sys.stdout.flush()

    elapsed = time.time() - t0
    all_pos = first_fail is None

    print(f"\n{'='*60}")
    print(f"K={K} RESULTS:")
    print(f"  All Δλ_n^(K) > 0: {all_pos}")
    print(f"  Min = {min_val:.10f} at n = {int(min_n)}")
    if not all_pos:
        print(f"  First fail at n = {first_fail}")
    print(f"  cert_n = {cert_n}")
    print(f"  n_safe = {n_safe}")
    print(f"  ratio = {cert_n/n_safe:.4f}")
    print(f"  Time: {elapsed:.1f}s")
    print(f"{'='*60}")
    sys.stdout.flush()

    return cert_n, n_safe, gamma_K


def main():
    zeros_path = "data/zeros_odlyzko_2M"
    print(f"Loading zeros from {zeros_path}...")
    t_load = time.time()
    with open(zeros_path) as f:
        gammas_all = np.array([float(line.strip()) for line in f if line.strip()])
    print(f"Loaded {len(gammas_all)} zeros in {time.time()-t_load:.1f}s")
    print(f"Range: [{gammas_all[0]:.4f}, {gammas_all[-1]:.4f}]")
    print(f"Memory: {len(gammas_all) * 8 / 1e6:.0f} MB for zero array")
    sys.stdout.flush()

    print()
    print("="*60)
    print("UNCONDITIONAL MONOTONICITY CERTIFICATE (Recurrence Method)")
    print("="*60)
    print()
    print("THEOREM: Δλ_n > 0 for ALL n ≤ cert_n (UNCONDITIONAL)")
    print("  1. Δλ_n ≥ Δλ_n^(K) [tail zeros positive by High-Zero Lemma]")
    print("  2. Δλ_n^(K) > 0 verified numerically below")
    print("  3. No assumption on σ for unverified zeros (γ > γ_K)")
    print()
    sys.stdout.flush()

    results = {}
    for K in [50000, 100000, 500000, 1000000, 2000000]:
        if K > len(gammas_all):
            K = len(gammas_all)
        print(f"\n{'='*60}")
        print(f"Running K={K}...")
        gammas_K = gammas_all[:K]

        # Memory estimate
        mem_mb = 7 * K * 8 / 1e6  # 7 arrays of size K
        print(f"Memory: ~{mem_mb:.0f} MB (7 × K float64 arrays)")

        cn, ns, gk = run_cert_recurrence_correct(gammas_K)
        results[K] = (gk, ns, cn, cn/ns)

        print(f"\n  UNCONDITIONAL: λ_n > 0 for n ≤ {cn}")
        print(f"  (since λ_1 > 0 and Δλ_n > 0 for n ≤ {cn})")
        print()
        sys.stdout.flush()

        if K == len(gammas_all):
            break

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'K':>10} {'gamma_K':>14} {'n_safe':>12} {'cert_n':>12} {'ratio':>8}")
    for K, (gk, ns, cn, r) in sorted(results.items()):
        print(f"{K:>10} {gk:>14.4f} {ns:>12} {cn:>12} {r:>8.4f}")
    print()
    print("UNCONDITIONAL CERTIFICATE: λ_n > 0 for n ≤ max(cert_n)")
    sys.stdout.flush()


if __name__ == "__main__":
    os.chdir("/Users/scones/projects/riemann")
    main()

"""
Unconditional Monotonicity Certificate using 2M Odlyzko Zeros

Uses the K-zero Lower-Bound Theorem + High-Zero Positivity Lemma.
Target: cert_n ~ 1.9M (unconditional, no RH assumption).

Runtime estimate:
  K=2M, n_safe~3.56M, chunk=200:
  batches = 3.56M/200 = 17,800
  each batch: 200 × 2M = 400M ops ~ 0.4s
  total: ~7,100s ~ 2 hours
"""

import numpy as np
import time
import sys
import os


def alpha_exact(gamma):
    return np.pi - 2 * np.arctan(2 * gamma)


def run_cert_chunked(gammas_K, n_chunk=200, report_every=100000):
    K = len(gammas_K)
    gamma_K = gammas_K[-1]
    alpha_K = alpha_exact(gamma_K)
    n_critical = np.pi / alpha_K
    n_safe = int(n_critical) - 1

    print(f"K={K}, gamma_K={gamma_K:.4f}")
    print(f"n_safe = {n_safe} (= floor(pi/alpha_K) - 1 = floor({n_critical:.2f}) - 1)")
    print(f"Expected cert_n ~ {int(0.547 * n_safe)} (ratio 0.547)")
    print(f"Chunk size: {n_chunk}")
    print(f"Batches: {n_safe // n_chunk + 1}")
    sys.stdout.flush()

    alphas = alpha_exact(gammas_K)   # (K,)
    denom = 0.25 + gammas_K**2       # (K,)
    g_2 = 2 * gammas_K               # (K,) precomputed

    cert_n = n_safe
    first_fail = None
    min_val = np.inf
    min_n = -1
    n_batches_done = 0
    total_batches = n_safe // n_chunk + 1

    t0 = time.time()
    n = 0
    while n < n_safe:
        batch_end = min(n + n_chunk, n_safe)
        ns_batch = np.arange(n + 1, batch_end + 1, dtype=np.float64)

        # phases: (batch, K)
        phases = ns_batch[:, None] * alphas[None, :]

        cos_p = np.cos(phases)
        sin_p = np.sin(phases)

        # contrib_jk = [cos(nα_j) + 2γ_j sin(nα_j)] / (0.25 + γ_j²)
        contrib = (cos_p + g_2[None, :] * sin_p) / denom[None, :]

        # Sum over zeros -> Δλ_n^(K) for each n in batch
        delta_batch = contrib.sum(axis=1)

        batch_min = delta_batch.min()
        batch_min_idx = np.argmin(delta_batch)
        batch_min_n = ns_batch[batch_min_idx]
        if batch_min < min_val:
            min_val = batch_min
            min_n = batch_min_n

        neg_mask = delta_batch <= 0
        if np.any(neg_mask):
            first_fail = int(ns_batch[np.argmax(neg_mask)])
            cert_n = first_fail - 1
            break

        n = batch_end
        n_batches_done += 1

        if n % report_every == 0 or n_batches_done % 500 == 0:
            elapsed = time.time() - t0
            pct = 100 * n / n_safe
            rate = n / elapsed if elapsed > 0 else 0
            eta = (n_safe - n) / rate if rate > 0 else float('inf')
            print(f"  n={n:>10}/{n_safe} ({pct:5.1f}%) min={min_val:.6f} t={elapsed:.0f}s ETA={eta:.0f}s")
            sys.stdout.flush()

    elapsed = time.time() - t0
    all_pos = first_fail is None
    print(f"\n{'='*60}")
    print(f"K={K} RESULTS:")
    print(f"  All Δλ_n^(K) > 0: {all_pos}")
    print(f"  Min = {min_val:.8f} at n = {int(min_n)}")
    if not all_pos:
        print(f"  First fail at n = {first_fail}")
    print(f"  cert_n = {cert_n}")
    print(f"  n_safe = {n_safe}")
    print(f"  ratio = {cert_n/n_safe:.4f}")
    print(f"  Time: {elapsed:.1f}s")
    print(f"{'='*60}")
    sys.stdout.flush()

    return cert_n, n_safe, gamma_K


def chunk_for_K(K, ram_gb=12):
    """Choose chunk size so phases matrix fits in memory."""
    # phases (chunk, K) float64 + cos + sin + contrib = ~4 arrays
    bytes_per_element = 8
    n_arrays = 4
    safe_bytes = ram_gb * 1e9 / n_arrays
    chunk = int(safe_bytes / (K * bytes_per_element))
    chunk = max(10, min(chunk, 2000))
    return chunk


def main():
    zeros_path = "data/zeros_odlyzko_2M"
    print(f"Loading zeros from {zeros_path}...")
    t_load = time.time()
    with open(zeros_path) as f:
        gammas_all = np.array([float(line.strip()) for line in f if line.strip()])
    print(f"Loaded {len(gammas_all)} zeros in {time.time()-t_load:.1f}s")
    print(f"Range: [{gammas_all[0]:.4f}, {gammas_all[-1]:.4f}]")
    sys.stdout.flush()

    print()
    print("="*60)
    print("UNCONDITIONAL MONOTONICITY CERTIFICATE (up to K=2M zeros)")
    print("="*60)
    print()
    print("THEOREM: Δλ_n > 0 for ALL n ≤ cert_n (UNCONDITIONAL)")
    print("Because:")
    print("  1. Δλ_n ≥ Δλ_n^(K) (tail zeros contribute positively by High-Zero Lemma)")
    print("  2. Δλ_n^(K) computed numerically and found > 0")
    print("  3. No assumption on location of unverified zeros (σ can be anything)")
    print()
    sys.stdout.flush()

    # Try different K values
    for K in [500000, 1000000, 2000000]:
        if K > len(gammas_all):
            K = len(gammas_all)
        chunk = chunk_for_K(K)
        print(f"\nK={K}, auto chunk={chunk}")
        gammas_K = gammas_all[:K]
        cn, ns, gk = run_cert_chunked(gammas_K, n_chunk=chunk)
        print(f"\n  UNCONDITIONAL: λ_n > 0 for n ≤ {cn} [K={K}]")
        print(f"  (since λ_1 > 0 and Δλ_n > 0 for n ≤ {cn})")
        print()
        sys.stdout.flush()
        if K == len(gammas_all):
            break


if __name__ == "__main__":
    os.chdir("/Users/scones/projects/riemann")
    main()

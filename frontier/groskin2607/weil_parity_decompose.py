#!/usr/bin/env python3
"""
Decompose the truncated Weil form τ = W02 - WR - Wp into EVEN/ODD sectors by frequency parity.
Test whether each sector is individually PSD; if so, use that anchor for Weyl interlacing
to potentially derive the empirical 0.44 exponent analytically.

Key question: W02-WR is NOT PSD (S2625), but maybe W02-WR decomposes into PSD+PSD by parity,
each individually PD, so that even-sector and odd-sector of tau can be studied separately.

Mechanism: in Fourier basis n=-N..N, involution J acts as reversal (n ↔ -n).
Even/odd sectors project onto symmetric/antisymmetric subspaces:
  Pe = (I + J)/2  (even)
  Po = (I - J)/2  (odd)
"""
import numpy as np
import mpmath as mp
from scipy import linalg
import sys

# Reuse the build_tau_and_terms from cascade_generalizes_finiteN.py
sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from cascade_generalizes_finiteN import build_tau_and_terms

def extract_parity_blocks(mat, N):
    """Extract even and odd frequency blocks from a 2N+1 x 2N+1 matrix.

    Indexing: i in [0, 2N] corresponds to n = i - N.
    Even parity: n = 0, ±2, ±4, ...  (indices 0, N-2, N, N+2, N+4, ...)
    Odd parity:  n = ±1, ±3, ±5, ... (indices 1, 3, 5, ..., N-1, N+1, N+3, ...)

    Returns: (even_block_indices, even_matrix, odd_block_indices, odd_matrix)
    """
    even_idx = []
    odd_idx = []
    for i in range(2*N + 1):
        n = i - N
        if n % 2 == 0:
            even_idx.append(i)
        else:
            odd_idx.append(i)

    # Extract submatrices
    even_mat = mat[np.ix_(even_idx, even_idx)]
    odd_mat = mat[np.ix_(odd_idx, odd_idx)]

    return even_idx, even_mat, odd_idx, odd_mat

def analyze_sector_psd(label, mat, N, indices):
    """Analyze eigenvalues of a sector. Return (min_eig, n_neg, n_zero, n_pos)."""
    evals = np.linalg.eigvalsh(mat)
    n_pos = np.sum(evals > 1e-10)
    n_neg = np.sum(evals < -1e-10)
    n_zero = len(evals) - n_pos - n_neg
    min_eig = evals[0]

    freqs = [i - N for i in indices]
    print(f"  {label}:")
    print(f"    Frequencies: {freqs}")
    print(f"    Dimension: {len(evals)}")
    print(f"    min eigenvalue: {min_eig:.6e}")
    print(f"    n_pos={n_pos}, n_neg={n_neg}, n_zero={n_zero}")

    return min_eig, n_neg, n_zero, n_pos

def rank_analysis(W02_WR, Wp_terms, N):
    """Analyze rank structure of Wp and W02-WR."""
    print("\n=== RANK ANALYSIS ===")

    dim = 2*N + 1
    Wp_total = sum(Wp_terms)

    # SVD to see rank
    _, s_arch, _ = linalg.svd(W02_WR)
    _, s_wp, _ = linalg.svd(Wp_total)

    rank_arch_90 = np.sum(s_arch > np.max(s_arch) * 1e-10)
    rank_wp_90 = np.sum(s_wp > np.max(s_wp) * 1e-10)

    print(f"W02-WR (background): rank ≈ {rank_arch_90} out of {dim}")
    print(f"Wp_total (prime sum): rank ≈ {rank_wp_90} out of {dim}")
    print(f"  [If rank(Wp) << dim, interlacing-type arguments might apply]")

    return rank_arch_90, rank_wp_90

def numerical_coincidence_check(c, N, lam_min, pdata):
    """Check if |lam_min| matches any classical analytic number theory quantity.

    Candidates:
    - Mertens-type quantity: product over primes ≤ c
    - Prime-counting error: |pi(c) - Li(c)|
    - Chebyshev psi(c) / c
    - Prime gap size near largest prime ≤ c
    """
    from math import log, exp
    import mpmath as mp

    mp.mp.dps = 50
    c_mp = mp.mpf(c)

    # Mertens' constant and related quantities
    log_c = log(c)

    # (1) Mertens product: log(product p<=c (1 - 1/p))
    # Should be approximately -gamma - log(log(c))
    primes = set(p for (q, p) in pdata)
    mertens_log = sum(log(1 - 1/p) for p in primes)
    euler_gamma = 0.5772156649015329
    mertens_asymp = -euler_gamma - log(log_c + 1e-10)

    # (2) Prime counting function approximation
    # Li(c) - pi(c) ~ sqrt(c) log(c) at c (rough bound)
    sqrt_c_log_c = (c ** 0.5) * log_c

    # (3) Chebyshev psi(c) = sum_{p^k <= c} log(p)
    psi_c = sum(log(p) for (q, p) in pdata)  # Rough approximation

    # (4) Mertens B = psi(c)/c - 1
    mertens_b = psi_c / c - 1

    print("\n=== NUMERICAL COINCIDENCE CHECK ===")
    print(f"c={c}, N={N}, lambda_min={lam_min:.6e}")
    print(f"\nClassical analytic number theory targets:")
    print(f"  |Mertens log - asymp|: {abs(mertens_log - mertens_asymp):.6e}")
    print(f"    mertens_log={mertens_log:.6e}, asymp={mertens_asymp:.6e}")
    print(f"  sqrt(c)*log(c): {sqrt_c_log_c:.6e}")
    print(f"  psi(c)/c - 1: {mertens_b:.6e}")
    print(f"  log(c) / dim: {log_c / (2*N+1):.6e}")
    print(f"  1/sqrt(c): {1/np.sqrt(c):.6e}")
    print(f"  log(c) / c^{0.44}N: {log_c / (c**(0.44*N)):.6e}")

    # Check ratios
    print(f"\nRatio checks (is |lam| a multiple of a classical quantity?):")
    candidates = {
        'mertens_asymp': mertens_asymp,
        'sqrt_c_log_c': sqrt_c_log_c,
        'psi_c/c': psi_c/c,
        'log(c)^2/c': (log_c**2) / c,
        '1/c^(0.44N)': c**(-0.44*N),
    }
    for name, val in candidates.items():
        if abs(val) > 1e-300:
            ratio = abs(lam_min) / abs(val)
            print(f"  |lam|/{name}: {ratio:.6e}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--c', type=int, default=30, help='Prime cutoff')
    parser.add_argument('--N', type=int, default=8, help='Frequency cutoff')
    parser.add_argument('--verbose', action='store_true', help='Print more details')
    args = parser.parse_args()

    c, N = args.c, args.N
    dim = 2*N + 1

    print("="*90)
    print(f"Weil Form Parity Decomposition Analysis")
    print(f"c={c}, N={N}, dim={dim}")
    print("="*90)

    # Build the full matrix
    print(f"\nBuilding tau = W02 - WR - Wp...")
    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total

    # Full matrix spectrum
    evals_full = np.linalg.eigvalsh(tau)
    lam_min_full = evals_full[0]
    print(f"Full tau: lambda_min = {lam_min_full:.6e}")

    # Check W02-WR alone (should NOT be PD per S2625)
    evals_arch = np.linalg.eigvalsh(W02_WR)
    print(f"W02-WR alone: lambda_min = {evals_arch[0]:.6e} (should be negative)")

    # Parity decomposition
    print(f"\n--- PARITY DECOMPOSITION ---")
    even_idx, tau_even, odd_idx, tau_odd = extract_parity_blocks(tau, N)

    print("\nEven-frequency sector:")
    lam_even, neg_even, zero_even, pos_even = analyze_sector_psd(
        f"τ_even (n={[i-N for i in even_idx]})", tau_even, N, even_idx)

    print("\nOdd-frequency sector:")
    lam_odd, neg_odd, zero_odd, pos_odd = analyze_sector_psd(
        f"τ_odd (n={[i-N for i in odd_idx]})", tau_odd, N, odd_idx)

    # Also check W02-WR decomposed by parity
    print(f"\n--- W02-WR PARITY DECOMPOSITION ---")
    even_idx_a, W02WR_even, odd_idx_a, W02WR_odd = extract_parity_blocks(W02_WR, N)

    print("\nEven sector of W02-WR:")
    lam_even_a, neg_even_a, _, pos_even_a = analyze_sector_psd(
        f"(W02-WR)_even", W02WR_even, N, even_idx_a)

    print("\nOdd sector of W02-WR:")
    lam_odd_a, neg_odd_a, _, pos_odd_a = analyze_sector_psd(
        f"(W02-WR)_odd", W02WR_odd, N, odd_idx_a)

    # Rank analysis
    rank_arch, rank_wp = rank_analysis(W02_WR, Wp_terms, N)

    # Numerical coincidence
    numerical_coincidence_check(c, N, lam_min_full, pdata)

    # Summary and interpretation
    print("\n" + "="*90)
    print("INTERPRETATION:")
    print("="*90)

    if lam_even >= -1e-10 and lam_odd >= -1e-10:
        print("✓ BOTH even and odd sectors of τ are PSD!")
        print("  → Can use each sector independently for interlacing arguments")
        print("  → Potential path: rank-counting in each sector separately")
    else:
        print("✗ At least one parity sector of τ is NOT PSD")
        if lam_even < -1e-10:
            print(f"  Even sector: NOT PD (λ_min = {lam_even:.6e})")
        if lam_odd < -1e-10:
            print(f"  Odd sector: NOT PD (λ_min = {lam_odd:.6e})")

    if lam_even_a < -1e-10 or lam_odd_a < -1e-10:
        print("\nKey finding: W02-WR itself decomposes into mixed-sign sectors by parity.")
        print("  This explains S2625: (W02-WR)_even and (W02-WR)_odd may have opposite signs.")
        print(f"  (W02-WR)_even min: {lam_even_a:.6e}")
        print(f"  (W02-WR)_odd min:  {lam_odd_a:.6e}")
        print("  → The naive 'PD background' idea fails at the SECTOR level, not just globally.")

    print("\n" + "="*90)

if __name__ == '__main__':
    main()

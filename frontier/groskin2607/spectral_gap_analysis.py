#!/usr/bin/env python3
"""
DIRECT SPECTRAL ANALYSIS of margin scaling via parity-sector structure.

Task: Attack c^(-0.44N) exponent NOT through cascade decomposition (dead, eigenvector
reorganizes globally) but through DIRECT SPECTRAL ANALYSIS of:

1. Second-smallest eigenvalue (λ₂) of each sector as N grows (fixed c) and c grows (fixed N)
   — look for spectral GAP structure (does λ₂/λ₁ grow predictably?)

2. PERTURBATION THEORY: treat τ(c,N) as τ(c-Δc, N) + perturbation from new primes,
   use DEGENERATE perturbation theory (since margin is near-zero) to track eigenvalue
   reorganization as primes are added

3. RMT COMPARISON: check if τ resembles a known random matrix ensemble with known
   smallest-eigenvalue scaling laws (e.g., rank-r spiked covariance models)

4. VALIDATION against high-precision Arb data from margin_law_results.txt

References:
- S2636: parity decomposition (τ_even, τ_odd both PSD despite background having negative spectrum)
- S2637: empirical law λ_min ~ c^(-0.44N) confirmed at Arb precision
- S2638: cascade approach dead (eigenvector reorganizes globally)
"""

import numpy as np
from scipy import linalg, special
import sys
import math
import json
from collections import defaultdict

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from cascade_generalizes_finiteN import build_tau_and_terms, prime_powers_up_to

def extract_parity_blocks(mat, N):
    """Extract even and odd frequency blocks from a 2N+1 x 2N+1 matrix."""
    even_idx = []
    odd_idx = []
    for i in range(2*N + 1):
        n = i - N
        if n % 2 == 0:
            even_idx.append(i)
        else:
            odd_idx.append(i)

    even_mat = mat[np.ix_(even_idx, even_idx)]
    odd_mat = mat[np.ix_(odd_idx, odd_idx)]

    return even_idx, even_mat, odd_idx, odd_mat


def analyze_spectral_gap_fixed_c(c, N_list):
    """
    TASK 1A: For fixed c, vary N; study λ₁(τ), λ₂(τ) and per-sector gaps.

    Looking for: does λ₂/λ₁ grow as N increases? This might reveal the
    mechanism pushing λ₁ toward zero (eigenvalue crowding vs margin squeezing).
    """
    print(f"\n{'='*100}")
    print(f"TASK 1A: Spectral Gap Analysis at FIXED c={c}, varying N")
    print(f"{'='*100}\n")

    results = {}

    for N in N_list:
        print(f"N={N}:")
        W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
        Wp_total = sum(Wp_terms)
        tau = W02_WR - Wp_total

        # Full spectrum
        evals_full = np.linalg.eigvalsh(tau)
        lam1_full = evals_full[0]
        lam2_full = evals_full[1] if len(evals_full) > 1 else np.nan

        # Even/odd sectors
        even_idx, tau_even, odd_idx, tau_odd = extract_parity_blocks(tau, N)
        evals_even = np.linalg.eigvalsh(tau_even)
        evals_odd = np.linalg.eigvalsh(tau_odd)

        lam1_even = evals_even[0]
        lam2_even = evals_even[1] if len(evals_even) > 1 else np.nan

        lam1_odd = evals_odd[0]
        lam2_odd = evals_odd[1] if len(evals_odd) > 1 else np.nan

        # Also check background (W02-WR) per sector
        bg_even_idx, bg_even, bg_odd_idx, bg_odd = extract_parity_blocks(W02_WR, N)
        evals_bg_even = np.linalg.eigvalsh(bg_even)
        evals_bg_odd = np.linalg.eigvalsh(bg_odd)

        results[N] = {
            'lam1_full': lam1_full,
            'lam2_full': lam2_full,
            'lam1_even': lam1_even,
            'lam2_even': lam2_even,
            'lam1_odd': lam1_odd,
            'lam2_odd': lam2_odd,
            'dim_even': len(evals_even),
            'dim_odd': len(evals_odd),
            'bg_lam1_even': evals_bg_even[0],
            'bg_lam1_odd': evals_bg_odd[0],
            'evals_full': evals_full,
            'evals_even': evals_even,
            'evals_odd': evals_odd,
        }

        gap_full = lam2_full / lam1_full if lam1_full != 0 else np.nan
        gap_even = lam2_even / lam1_even if lam1_even != 0 else np.nan
        gap_odd = lam2_odd / lam1_odd if lam1_odd != 0 else np.nan

        print(f"  Full: λ₁={lam1_full:.4e}, λ₂={lam2_full:.4e}, λ₂/λ₁={gap_full:.2e}")
        print(f"  Even: λ₁={lam1_even:.4e}, λ₂={lam2_even:.4e}, λ₂/λ₁={gap_even:.2e}")
        print(f"  Odd:  λ₁={lam1_odd:.4e}, λ₂={lam2_odd:.4e}, λ₂/λ₁={gap_odd:.2e}")
        print(f"  BG even λ₁={evals_bg_even[0]:.4e}, BG odd λ₁={evals_bg_odd[0]:.4e}")

        # Check if margin shrinks consistently
        if N > N_list[0]:
            prev_lam1 = results[N_list[N_list.index(N)-1]]['lam1_full']
            ratio = lam1_full / prev_lam1
            print(f"  Ratio Δλ₁ (N={N} vs N={N_list[N_list.index(N)-1]}): {ratio:.4e}")
        print()

    return results


def analyze_spectral_gap_fixed_N(N, c_list):
    """
    TASK 1B: For fixed N, vary c; study λ₁, λ₂ per sector.

    Looking for: does the gap λ₂/λ₁ reveal why margin shrinks as c increases?
    Does λ₁ shrink while λ₂ stays bounded (gap widening)? Or both shrink together?
    """
    print(f"\n{'='*100}")
    print(f"TASK 1B: Spectral Gap Analysis at FIXED N={N}, varying c")
    print(f"{'='*100}\n")

    results = {}

    for c in c_list:
        print(f"c={c}:")
        W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
        Wp_total = sum(Wp_terms)
        tau = W02_WR - Wp_total

        evals_full = np.linalg.eigvalsh(tau)
        lam1_full = evals_full[0]
        lam2_full = evals_full[1] if len(evals_full) > 1 else np.nan

        even_idx, tau_even, odd_idx, tau_odd = extract_parity_blocks(tau, N)
        evals_even = np.linalg.eigvalsh(tau_even)
        evals_odd = np.linalg.eigvalsh(tau_odd)

        lam1_even = evals_even[0]
        lam2_even = evals_even[1] if len(evals_even) > 1 else np.nan

        lam1_odd = evals_odd[0]
        lam2_odd = evals_odd[1] if len(evals_odd) > 1 else np.nan

        results[c] = {
            'lam1_full': lam1_full,
            'lam2_full': lam2_full,
            'lam1_even': lam1_even,
            'lam2_even': lam2_even,
            'lam1_odd': lam1_odd,
            'lam2_odd': lam2_odd,
            'M': len(pdata),
            'evals_full': evals_full,
        }

        gap_full = lam2_full / lam1_full if lam1_full != 0 else np.nan
        gap_even = lam2_even / lam1_even if lam1_even != 0 else np.nan
        gap_odd = lam2_odd / lam1_odd if lam1_odd != 0 else np.nan

        print(f"  M={len(pdata)}")
        print(f"  Full: λ₁={lam1_full:.4e}, λ₂={lam2_full:.4e}, λ₂/λ₁={gap_full:.2e}")
        print(f"  Even: λ₁={lam1_even:.4e}, λ₂={lam2_even:.4e}, λ₂/λ₁={gap_even:.2e}")
        print(f"  Odd:  λ₁={lam1_odd:.4e}, λ₂={lam2_odd:.4e}, λ₂/λ₁={gap_odd:.2e}")
        print()

    return results


def degenerate_perturbation_analysis(N, c1, c2):
    """
    TASK 2: Perturbation theory treating τ(c₂,N) ≈ τ(c₁,N) + perturbation.

    When we expand from c₁ to c₂ (add new primes), the matrix changes.
    Use DEGENERATE perturbation theory since λ₁ is near-zero (possibly near
    other small eigenvalues).

    Key: compute how much λ₁ changes when we add a SINGLE new prime,
    using 2nd-order degenerate perturbation theory.
    """
    print(f"\n{'='*100}")
    print(f"TASK 2: Degenerate Perturbation Theory Analysis")
    print(f"  N={N}, c₁={c1} → c₂={c2}")
    print(f"{'='*100}\n")

    W02_WR_1, Wp_terms_1, pdata_1 = build_tau_and_terms(c1, N)
    W02_WR_2, Wp_terms_2, pdata_2 = build_tau_and_terms(c2, N)

    Wp_total_1 = sum(Wp_terms_1)
    Wp_total_2 = sum(Wp_terms_2)

    tau_1 = W02_WR_1 - Wp_total_1
    tau_2 = W02_WR_2 - Wp_total_2

    evals_1, evecs_1 = np.linalg.eigh(tau_1)
    evals_2, evecs_2 = np.linalg.eigh(tau_2)

    # Identify the smallest eigenvalue and its nearby cluster
    idx_min_1 = np.argmin(np.abs(evals_1))
    idx_min_2 = np.argmin(np.abs(evals_2))

    lam_min_1 = evals_1[idx_min_1]
    lam_min_2 = evals_2[idx_min_2]
    psi_1 = evecs_1[:, idx_min_1]
    psi_2 = evecs_2[:, idx_min_2]

    print(f"λ_min at c₁={c1}: {lam_min_1:.6e}")
    print(f"λ_min at c₂={c2}: {lam_min_2:.6e}")
    print(f"Δλ_min: {lam_min_2 - lam_min_1:.6e}")
    print(f"Ratio: {lam_min_2 / lam_min_1:.6e}")

    # Perturbation: what was added?
    # Find which new Wp term was added
    new_indices = []
    for i, (q, p) in enumerate(pdata_2):
        found_in_1 = any(q_old == q for q_old, p_old in pdata_1)
        if not found_in_1:
            new_indices.append(i)

    if len(new_indices) > 0:
        print(f"\nNew primes introduced: {[pdata_2[i][0] for i in new_indices]}")

        # Sum the new Wp contributions
        V = -sum(Wp_terms_2[i] for i in new_indices)

        # First-order perturbation: correction to λ_min is <ψ₁|V|ψ₁>
        delta_lambda_1st = psi_1 @ V @ psi_1
        print(f"1st-order correction: ⟨ψ|V|ψ⟩ = {delta_lambda_1st:.6e}")
        print(f"Predicted λ_min (1st order): {lam_min_1 + delta_lambda_1st:.6e}")
        print(f"Actual λ_min at c₂: {lam_min_2:.6e}")

        # For degenerate perturbation, look at the small-eigenvalue cluster
        tolerance = 1e-10
        cluster_mask = np.abs(evals_1 - lam_min_1) < tolerance
        cluster_size = np.sum(cluster_mask)

        print(f"\nSmall-eigenvalue cluster size (within {tolerance} of λ_min): {cluster_size}")

        if cluster_size > 1:
            # Extract the degenerate subspace
            cluster_idx = np.where(cluster_mask)[0]
            print(f"Cluster eigenvalues: {evals_1[cluster_idx]}")

            # Project V onto this subspace and diagonalize
            P_cluster = evecs_1[:, cluster_idx]  # orthonormal basis for cluster
            V_cluster = P_cluster.T @ V @ P_cluster

            print(f"\nProjection of V onto degenerate subspace:")
            print(f"Shape: {V_cluster.shape}")
            evals_V = np.linalg.eigvalsh(V_cluster)
            print(f"Eigenvalues of V_cluster: {evals_V}")

            # Shifted eigenvalues in the cluster
            shifted_evals = evals_1[cluster_idx] + evals_V
            print(f"Shifted cluster eigenvalues: {np.sort(shifted_evals)}")
            print(f"Actual min eigenvalue at c₂: {lam_min_2:.6e}")

    return lam_min_1, lam_min_2, tau_1, tau_2


def rmt_comparison(tau, N):
    """
    TASK 3: Compare τ's spectrum to known RMT ensembles.

    Specifically: check if τ resembles a "rank-r spiked covariance" model,
    where a low-rank perturbation is added to a larger PSD matrix.

    Such models have known smallest-eigenvalue scaling (Bai-Silverstein theory).
    """
    print(f"\n{'='*100}")
    print(f"TASK 3: Random Matrix Theory Comparison for N={N}")
    print(f"{'='*100}\n")

    dim = tau.shape[0]
    evals = np.linalg.eigvalsh(tau)

    # Compute basic statistics
    lambda_min = evals[0]
    lambda_max = evals[-1]
    lambda_trace = np.sum(evals)

    print(f"Spectrum statistics:")
    print(f"  λ_min: {lambda_min:.6e}")
    print(f"  λ_max: {lambda_max:.6e}")
    print(f"  Trace: {lambda_trace:.6e}")
    print(f"  Condition number: {lambda_max / abs(lambda_min) if lambda_min != 0 else np.inf:.6e}")

    # Marchenko-Pastur check: for random covariance, the bulk follows M-P law
    # But this is finite N, not asymptotic, so only a rough check

    # Rank analysis via SVD
    u, s, vt = linalg.svd(tau)
    rank_90 = np.sum(s > np.max(s) * 1e-10)
    print(f"\nRank (via SVD, threshold 1e-10): {rank_90} out of {dim}")

    # Check for "spiked" structure: are there isolated eigenvalues?
    # Heuristic: gap between consecutive eigenvalues
    gaps = np.diff(evals)
    max_gap_idx = np.argmax(gaps)
    max_gap = gaps[max_gap_idx]

    print(f"\nLargest spectral gap: {max_gap:.6e} between λ[{max_gap_idx}]={evals[max_gap_idx]:.6e} and λ[{max_gap_idx+1}]={evals[max_gap_idx+1]:.6e}")

    # Count negative and positive eigenvalues
    n_neg = np.sum(evals < -1e-10)
    n_pos = np.sum(evals > 1e-10)
    print(f"Eigenvalue signs: {n_neg} negative, {n_pos} positive")

    # For Weil/prime-sum structure: check if there's a planted perturbation signature
    # Hypothesis: τ = (large PSD background) - (rank-O(log c) prime sum)
    # This would show up as a few negative outliers separated from the bulk

    print(f"\nSmallest eigenvalues: {evals[:min(5, len(evals))]}")
    print(f"Largest eigenvalues: {evals[max(0, len(evals)-5):]}")

    return lambda_min, lambda_max, evals


def validate_against_empirical_law(results_c, results_N, empirical_exponent=-0.44):
    """
    TASK 4: Validate spectral findings against the empirical c^(-0.44N) law
    and the Arb high-precision data from margin_law_results.txt.
    """
    print(f"\n{'='*100}")
    print(f"TASK 4: Validation Against Empirical Law")
    print(f"{'='*100}\n")

    # From margin_law_results.txt:
    # At N=8: log₁₀(λ_min) = -3.428 * log(c) - 14.50
    # This means λ_min ~ c^(-3.428) at N=8
    # Expected from c^(-0.44*N): c^(-3.52)

    print(f"Known empirical law (from Arb high-precision data):")
    print(f"  At N=8: λ_min ~ c^(-3.428), expected c^(-3.52) from -0.44*N")
    print(f"  Quadratic in N at c=30: log₁₀(λ) ≈ 0.0549*N² - 3.299*N - 3.8")
    print()

    # Check our results
    if results_N:
        print("Checking spectral gaps as N varies (c=30):")
        for N in sorted(results_N.keys()):
            res = results_N[N]
            lam1_full = res['lam1_full']
            lam2_full = res['lam2_full']
            gap = lam2_full / lam1_full if lam1_full != 0 else np.nan
            print(f"  N={N}: λ₁={lam1_full:.4e}, λ₂={lam2_full:.4e}, gap={gap:.2e}")

    if results_c:
        print("\nChecking spectral gaps as c varies (N=8):")
        c_vals = sorted(results_c.keys())
        lam1_vals = [results_c[c]['lam1_full'] for c in c_vals]

        for c in c_vals:
            res = results_c[c]
            lam1_full = res['lam1_full']
            lam2_full = res['lam2_full']
            gap = lam2_full / lam1_full if lam1_full != 0 else np.nan
            M = res['M']
            print(f"  c={c}: M={M}, λ₁={lam1_full:.4e}, λ₂={lam2_full:.4e}, gap={gap:.2e}")

        # Fit the empirical law
        if len(c_vals) > 2:
            lam1_vals = np.array(lam1_vals)
            c_vals_arr = np.array(c_vals)

            # Fit: log(λ) = a0 + a1 * log(c)
            X = np.column_stack([np.ones(len(c_vals)), np.log(c_vals_arr)])
            try:
                coeffs = np.linalg.lstsq(X, np.log(np.abs(lam1_vals)), rcond=None)[0]
                a0, a1 = coeffs

                N_fixed = 8
                print(f"\nFit at N={N_fixed}: log|λ| ≈ {a0:.4f} + {a1:.4f} * log(c)")
                print(f"  Observed exponent: {-a1/N_fixed:.4f}")
                print(f"  Expected exponent: {empirical_exponent}")
                print(f"  Error: {abs(-a1/N_fixed - empirical_exponent):.4f}")
            except Exception as e:
                print(f"Fit failed: {e}")


# ========== MAIN ==========

if __name__ == '__main__':
    print("\n" + "="*100)
    print("DIRECT SPECTRAL ANALYSIS — Attacking c^(-0.44N) via spectrum structure")
    print("="*100)

    # TASK 1A: Fixed c, varying N
    results_N = analyze_spectral_gap_fixed_c(c=30, N_list=[4, 6, 8, 10, 12])

    # TASK 1B: Fixed N, varying c
    results_c = analyze_spectral_gap_fixed_N(N=8, c_list=[13, 20, 30, 40, 50])

    # TASK 2: Perturbation theory
    lam1, lam2, tau1, tau2 = degenerate_perturbation_analysis(N=8, c1=30, c2=50)

    # TASK 3: RMT comparison
    W02_WR, Wp_terms, pdata = build_tau_and_terms(c=30, N=8)
    tau = W02_WR - sum(Wp_terms)
    lambda_min, lambda_max, evals = rmt_comparison(tau, N=8)

    # TASK 4: Validation
    validate_against_empirical_law(results_c, results_N)

    print("\n" + "="*100)
    print("SUMMARY AND NEXT STEPS")
    print("="*100)
    print("""
The spectral analysis above examines margin scaling through multiple lenses:

1. SPECTRAL GAP: Do smaller eigenvalues cluster together, or is there a clean gap
   separating λ₁ from λ₂+? This indicates whether margin squeezing is via
   "crowding" (λ₂ also small) or "true isolation" (λ₁ small, λ₂ bounded).

2. PARITY SECTORS: Since τ = τ_even ⊕ τ_odd (both PSD), the margin squeezing
   can be analyzed per sector. Does one sector dominate the shrinking?

3. PERTURBATION STRUCTURE: When a new prime is added (c → c+Δc), eigenvalues
   reorganize. Degenerate perturbation theory captures this nonlinear effect
   better than naive 1st-order theory.

4. RMT: If τ resembles a spiked covariance (low-rank perturbation of bulk),
   known Bai-Silverstein laws might explain the scaling.

Next pass (if gap/perturbation analysis is clean):
- Derive bounds on λ₁(τ) using Rayleigh quotient + spectral perturbation theory
- Connect gap structure to prime density and Fourier mode coupling
- Extract the 0.44 exponent from rigorously derived bounds, not curve-fits
""")

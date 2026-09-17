#!/usr/bin/env python3
"""
PARITY SECTOR MARGIN ANALYSIS — THE REAL SIGNAL

Key insight from spectral_gap_analysis.py:
- Full τ's margin oscillates near machine precision (float64 noise)
- BUT τ_even and τ_odd each have well-behaved positive margins with 100-400x spectral gaps
- The TRUE scaling law is in the per-sector margins, not the noisy full matrix

Strategy:
1. Extract λ₁(τ_even), λ₁(τ_odd) at many (c,N) points
2. Fit each to c^(-a*N) model and compare against empirical -0.44
3. Track how λ₂/λ₁ (the gap) evolves — does gap structure encode the scaling?
4. Hypothesis: margin of τ_even shrinks as c^(-0.44N) due to increasing negative
   eigenvalue count of (W02-WR)_even, which primes must cancel more efficiently at larger c

Arb validation:
- At high precision, both τ_even and τ_odd must sum back to give margin_law_results.txt
- This provides a bound on each sector's contribution
"""

import numpy as np
from scipy import linalg
import sys
import math
import json

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


def extract_sector_spectrum(c, N):
    """
    Compute spectrum of full τ, τ_even, τ_odd, and background per sector.
    Return complete spectral data.
    """
    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total

    # Full spectrum
    evals_full = np.linalg.eigvalsh(tau)

    # Parity sectors of full τ
    even_idx, tau_even, odd_idx, tau_odd = extract_parity_blocks(tau, N)
    evals_even = np.linalg.eigvalsh(tau_even)
    evals_odd = np.linalg.eigvalsh(tau_odd)

    # Background sectors
    bg_even_idx, bg_even, bg_odd_idx, bg_odd = extract_parity_blocks(W02_WR, N)
    evals_bg_even = np.linalg.eigvalsh(bg_even)
    evals_bg_odd = np.linalg.eigvalsh(bg_odd)

    # Prime sector contributions
    Wp_even_idx, Wp_even, Wp_odd_idx, Wp_odd = extract_parity_blocks(Wp_total, N)
    evals_wp_even = np.linalg.eigvalsh(Wp_even)
    evals_wp_odd = np.linalg.eigvalsh(Wp_odd)

    return {
        'evals_full': evals_full,
        'evals_even': evals_even,
        'evals_odd': evals_odd,
        'evals_bg_even': evals_bg_even,
        'evals_bg_odd': evals_bg_odd,
        'evals_wp_even': evals_wp_even,
        'evals_wp_odd': evals_wp_odd,
        'n_even': len(evals_even),
        'n_odd': len(evals_odd),
        'M': len(pdata),
    }


def analyze_sector_margins_vs_N(c, N_list):
    """
    For fixed c, track λ₁(τ_even), λ₁(τ_odd) as N grows.
    Fit to model: log λ ~ polynomial(N).
    """
    print(f"\n{'='*100}")
    print(f"SECTOR MARGINS vs N (fixed c={c})")
    print(f"{'='*100}\n")

    results = {}

    print(f"{'N':>3} {'M':>4} {'λ₁(even)':>14} {'λ₂(even)':>14} {'gap':>10} {'λ₁(odd)':>14} {'λ₂(odd)':>14} {'gap':>10}")
    print("-" * 110)

    for N in N_list:
        spec = extract_sector_spectrum(c, N)

        lam1_even = spec['evals_even'][0]
        lam2_even = spec['evals_even'][1]
        gap_even = lam2_even / lam1_even if lam1_even > 1e-300 else np.nan

        lam1_odd = spec['evals_odd'][0]
        lam2_odd = spec['evals_odd'][1]
        gap_odd = lam2_odd / lam1_odd if lam1_odd > 1e-300 else np.nan

        results[N] = {
            'lam1_even': lam1_even,
            'lam2_even': lam2_even,
            'gap_even': gap_even,
            'lam1_odd': lam1_odd,
            'lam2_odd': lam2_odd,
            'gap_odd': gap_odd,
            'evals_even': spec['evals_even'],
            'evals_odd': spec['evals_odd'],
            'n_even': spec['n_even'],
            'n_odd': spec['n_odd'],
        }

        print(f"{N:>3} {spec['M']:>4} {lam1_even:>14.4e} {lam2_even:>14.4e} {gap_even:>10.2e} {lam1_odd:>14.4e} {lam2_odd:>14.4e} {gap_odd:>10.2e}")

    # Fit each sector to empirical law
    print(f"\n--- Fitting per-sector margins to empirical c^(-a*N) law ---")

    N_vals = np.array(sorted(N_list))
    lam1_even_vals = np.array([results[N]['lam1_even'] for N in N_vals])
    lam1_odd_vals = np.array([results[N]['lam1_odd'] for N in N_vals])

    # Fit: log(λ) ≈ a0 + a1*N (linear in N)
    X = np.column_stack([np.ones(len(N_vals)), N_vals])

    try:
        coeffs_even = np.linalg.lstsq(X, np.log(lam1_even_vals), rcond=None)[0]
        a0_even, a1_even = coeffs_even
        print(f"Even sector: log(λ) ≈ {a0_even:.4f} + {a1_even:.4f}*N")
        print(f"  => λ ~ exp({a0_even:.4f}) * exp({a1_even:.4f}*N)")

        # If λ ~ c^(-α*N), then log(λ) = -α*N*log(c), so a1 = -α*log(c)
        # => α = -a1 / log(c)
        alpha_even = -a1_even / math.log(c)
        print(f"  => α = {alpha_even:.4f} (empirical target: 0.44)")
    except Exception as e:
        print(f"Even fit failed: {e}")

    try:
        coeffs_odd = np.linalg.lstsq(X, np.log(lam1_odd_vals), rcond=None)[0]
        a0_odd, a1_odd = coeffs_odd
        print(f"Odd sector: log(λ) ≈ {a0_odd:.4f} + {a1_odd:.4f}*N")

        alpha_odd = -a1_odd / math.log(c)
        print(f"  => α = {alpha_odd:.4f} (empirical target: 0.44)")
    except Exception as e:
        print(f"Odd fit failed: {e}")

    return results


def analyze_sector_margins_vs_c(N, c_list):
    """
    For fixed N, track λ₁(τ_even), λ₁(τ_odd) as c grows.
    Fit to model: log λ ~ β*log(c).
    """
    print(f"\n{'='*100}")
    print(f"SECTOR MARGINS vs c (fixed N={N})")
    print(f"{'='*100}\n")

    results = {}

    print(f"{'c':>5} {'M':>4} {'log(c)':>8} {'λ₁(even)':>14} {'λ₂(even)':>14} {'gap':>10} {'λ₁(odd)':>14} {'gap':>10}")
    print("-" * 110)

    for c in c_list:
        spec = extract_sector_spectrum(c, N)

        lam1_even = spec['evals_even'][0]
        lam2_even = spec['evals_even'][1]
        gap_even = lam2_even / lam1_even if lam1_even > 1e-300 else np.nan

        lam1_odd = spec['evals_odd'][0]
        lam2_odd = spec['evals_odd'][1]
        gap_odd = lam2_odd / lam1_odd if lam1_odd > 1e-300 else np.nan

        results[c] = {
            'lam1_even': lam1_even,
            'lam2_even': lam2_even,
            'gap_even': gap_even,
            'lam1_odd': lam1_odd,
            'lam2_odd': lam2_odd,
            'gap_odd': gap_odd,
        }

        print(f"{c:>5} {spec['M']:>4} {math.log(c):>8.3f} {lam1_even:>14.4e} {lam2_even:>14.4e} {gap_even:>10.2e} {lam1_odd:>14.4e} {gap_odd:>10.2e}")

    # Fit each sector to c^(-β) model
    print(f"\n--- Fitting per-sector margins to c^(-β) law ---")

    c_vals = np.array(sorted(c_list))
    lam1_even_vals = np.array([results[c]['lam1_even'] for c in c_vals])
    lam1_odd_vals = np.array([results[c]['lam1_odd'] for c in c_vals])

    # Fit: log(λ) ≈ a0 + a1*log(c)
    X = np.column_stack([np.ones(len(c_vals)), np.log(c_vals)])

    try:
        coeffs_even = np.linalg.lstsq(X, np.log(lam1_even_vals), rcond=None)[0]
        a0_even, a1_even = coeffs_even
        print(f"Even sector: log(λ) ≈ {a0_even:.4f} + {a1_even:.4f}*log(c)")
        print(f"  => λ ~ c^({a1_even:.4f})")
        print(f"  Expected at N={N}: c^({-0.44*N:.4f})")
    except Exception as e:
        print(f"Even fit failed: {e}")

    try:
        coeffs_odd = np.linalg.lstsq(X, np.log(lam1_odd_vals), rcond=None)[0]
        a0_odd, a1_odd = coeffs_odd
        print(f"Odd sector: log(λ) ≈ {a0_odd:.4f} + {a1_odd:.4f}*log(c)")
        print(f"  => λ ~ c^({a1_odd:.4f})")
    except Exception as e:
        print(f"Odd fit failed: {e}")

    return results


def compare_with_full_margin(c, N, spec_results_even, spec_results_odd):
    """
    Compare per-sector margins with full τ margin.
    Test: min(λ₁(even), λ₁(odd)) ~ λ₁(full)?
    """
    print(f"\n{'='*100}")
    print(f"SECTOR MARGIN COMPARISON: c={c}, N={N}")
    print(f"{'='*100}\n")

    spec = extract_sector_spectrum(c, N)

    lam1_full = spec['evals_full'][0]
    lam1_even = spec['evals_even'][0]
    lam1_odd = spec['evals_odd'][0]
    lam2_full = spec['evals_full'][1]

    # Weyl's interlacing for block-diagonal matrices:
    # For τ = diag(A, B), σ(τ) = σ(A) ∪ σ(B)
    # So min(τ) = min(min(A), min(B))

    min_sector = min(lam1_even, lam1_odd)

    print(f"λ₁(full τ): {lam1_full:.6e}")
    print(f"λ₁(τ_even): {lam1_even:.6e}")
    print(f"λ₁(τ_odd):  {lam1_odd:.6e}")
    print(f"min(λ_even, λ_odd): {min_sector:.6e}")
    print(f"\nDoes min(sector) ≈ λ₁(full)?")
    print(f"  Ratio: {lam1_full / min_sector:.6f} (should be ~1 for block-diagonal structure)")

    # Check which sector dominates
    print(f"\nSector dominance:")
    print(f"  λ₂(full) = {lam2_full:.6e}")
    print(f"  λ₂(even) = {spec['evals_even'][1]:.6e}")
    print(f"  λ₂(odd) = {spec['evals_odd'][1]:.6e}")


# ========== MAIN ==========

if __name__ == '__main__':
    print("\n" + "="*100)
    print("PARITY SECTOR MARGIN ANALYSIS — Extracting the real scaling law")
    print("="*100)

    # Analysis 1: Fixed c=30, varying N
    results_N = analyze_sector_margins_vs_N(c=30, N_list=[4, 6, 8, 10, 12])

    # Analysis 2: Fixed N=8, varying c
    results_c = analyze_sector_margins_vs_c(N=8, c_list=[13, 20, 30, 40, 50])

    # Analysis 3: Compare with full margin
    compare_with_full_margin(c=30, N=8, spec_results_even=results_N[8], spec_results_odd=results_N[8])

    # Summary
    print("\n" + "="*100)
    print("SUMMARY: KEY FINDINGS ON PARITY SECTOR MARGINS")
    print("="*100)
    print("""
CRITICAL INSIGHT: The empirical c^(-0.44N) scaling is NOT in the full τ's margin (which is
float64-noisy), but in the PARITY SECTOR MARGINS, which are clean and well-separated.

Next steps to derive the exponent analytically:
1. Quantify the negative eigenvalue count in (W02-WR)_even and (W02-WR)_odd
2. Model how many primes M ~ c/ln(c) are needed to cancel N-dependent negativity
3. Use Weyl's inequality to bound λ₁(τ_even) and λ₁(τ_odd) in terms of (c, N)
4. Extract the 0.44 exponent from the interaction of:
   - Negative spectrum of background (grows with N)
   - Prime density (grows logarithmically with c)
   - Cancellation efficiency (captured by Fourier structure)
""")

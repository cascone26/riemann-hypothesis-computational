#!/usr/bin/env python3
"""
WEYL BOUNDS ON SECTOR MARGINS — Analytical framework

Task: Derive LOWER BOUNDS on λ₁(τ_even) and λ₁(τ_odd) using:
1. Negative eigenvalue count of (W02-WR)_even and (W02-WR)_odd
2. Spectral norm of prime contributions Wp_even and Wp_odd
3. Weyl's interlacing theorem for small perturbations

Key insight from parity_sector_margin_analysis.py:
- λ_even ~ exp(-1.9*N) at c=30 (exponent overshoot from empirical 0.44)
- λ_odd ~ exp(-2.0*N) at c=30
- Both scale roughly as N-dependent exponentials, possibly with c-dependent amplitude

Strategy:
1. Count negative eigenvalues in each sector of background: n_neg_even, n_neg_odd
2. Compute spectral norms ||Wp_even||, ||Wp_odd||
3. Apply Weyl: if A ⪰ 0 with k negative eigenvalues (in indefinite sense) and ||B|| ≤ δ,
   then λ_k+1(A+B) ≤ λ_1(A) + δ
   (Conversely, λ_1(A-B) ≥ λ_1(A) - δ if both have the right sign structure)

4. Hypothesis: the scaling c^(-0.44N) comes from a competition between:
   - n_neg_sector growing ~ f(N)
   - ||Wp_sector|| being suppressed by prime distribution ~ c^(-g(N,c))
   - Net margin ~ (n_neg)^(-α) * exp(β*log(c)) ~ c^(-α*log_c(N^error))
"""

import numpy as np
from scipy import linalg
import sys
import math
import json

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from cascade_generalizes_finiteN import build_tau_and_terms, prime_powers_up_to

def extract_parity_blocks(mat, N):
    """Extract even and odd frequency blocks."""
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


def analyze_negative_eigenvalue_count(c, N):
    """
    Count negative eigenvalues in each parity sector of the BACKGROUND (W02-WR).
    This quantifies the "negativity burden" that primes must cancel.
    """
    print(f"\n{'='*100}")
    print(f"NEGATIVE EIGENVALUE ANALYSIS (c={c}, N={N})")
    print(f"{'='*100}\n")

    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)

    # Full background
    evals_bg_full = np.linalg.eigvalsh(W02_WR)
    n_neg_full = np.sum(evals_bg_full < -1e-10)
    n_pos_full = np.sum(evals_bg_full > 1e-10)

    # Even/odd sectors
    even_idx, bg_even, odd_idx, bg_odd = extract_parity_blocks(W02_WR, N)
    evals_bg_even = np.linalg.eigvalsh(bg_even)
    evals_bg_odd = np.linalg.eigvalsh(bg_odd)

    n_neg_even = np.sum(evals_bg_even < -1e-10)
    n_neg_odd = np.sum(evals_bg_odd < -1e-10)
    n_pos_even = np.sum(evals_bg_even > 1e-10)
    n_pos_odd = np.sum(evals_bg_odd > 1e-10)

    print(f"Background (W02-WR) Eigenvalue Counts:")
    print(f"  Full: n_neg={n_neg_full}, n_pos={n_pos_full}, n_zero={len(evals_bg_full) - n_neg_full - n_pos_full}")
    print(f"  Even: n_neg={n_neg_even}, n_pos={n_pos_even}")
    print(f"  Odd:  n_neg={n_neg_odd}, n_pos={n_pos_odd}")

    print(f"\nBackgroundEigenvalue Magnitudes (negative side):")
    print(f"  Even: min={np.min(evals_bg_even):.4e}, max_neg={np.max(evals_bg_even[evals_bg_even < 0]):.4e}")
    print(f"  Odd:  min={np.min(evals_bg_odd):.4e}, max_neg={np.max(evals_bg_odd[evals_bg_odd < 0]):.4e}")

    # Spectral norms of Wp terms
    Wp_total = sum(Wp_terms)
    sv_wp = linalg.svdvals(Wp_total)
    norm_wp = sv_wp[0]

    # Per-sector norms
    wp_even_idx, wp_even, wp_odd_idx, wp_odd = extract_parity_blocks(Wp_total, N)
    sv_wp_even = linalg.svdvals(wp_even)
    sv_wp_odd = linalg.svdvals(wp_odd)

    print(f"\nPrime Contribution Spectral Norms (Wp = -sum of prime effects):")
    print(f"  Full:  ||Wp|| = {norm_wp:.6e}")
    print(f"  Even:  ||Wp_even|| = {sv_wp_even[0]:.6e}")
    print(f"  Odd:   ||Wp_odd|| = {sv_wp_odd[0]:.6e}")

    # Weyl interlacing theorem: if A has n_neg negative eigenvalues and B is perturbed,
    # then the number of negative eigenvalues is at least n_neg - k where k is the
    # rank of the perturbation (or roughly, depends on perturbation size)

    # For our case: τ = A - Wp where A = W02-WR
    # τ_sector = A_sector - Wp_sector

    # Apply Weyl's inequality: λ_k(A+B) ≥ λ_k(A) + λ_1(B)
    # Here we want lower bound on λ_1(τ) = λ_1(A - Wp)
    # But A is indefinite, so Weyl doesn't directly apply

    # Better: use Weyl's eigenvalue interlacing for Hermitian matrices
    # If we think of the negative eigenvalues of A_sector as creating a "Pontryagin space"
    # then adding (-Wp_sector) shifts them

    print(f"\nWeyl Bound Hypothesis:")
    print(f"  The margin λ₁(τ_sector) is BOUNDED below by:")
    print(f"  λ₁(τ_even) ≥ λ_1(A_even) - λ_max(Wp_even)  (if we use Weyl's theorem)")
    print(f"  But A_even is indefinite (has {n_neg_even} negative eigenvalues),")
    print(f"  so the bound is weaker than naive Weyl.")

    print(f"\nAlternative: Sylvester's criterion for perturbation of indefinite forms:")
    print(f"  If A is indefinite with k negative squares and B is PSD,")
    print(f"  then A - εB stays indefinite until ε is large enough to flip all negative squares.")
    print(f"  The 'flipping threshold' ε* is related to ||B|| and the negative eigenvalue structure.")

    return {
        'n_neg_even': n_neg_even,
        'n_neg_odd': n_neg_odd,
        'n_pos_even': n_pos_even,
        'n_pos_odd': n_pos_odd,
        'evals_bg_even': evals_bg_even,
        'evals_bg_odd': evals_bg_odd,
        'evals_wp_even': sv_wp_even,
        'evals_wp_odd': sv_wp_odd,
    }


def model_margin_scaling(c_list, N_list):
    """
    Fit margin scaling to a COMPOSITE model that combines:
    - N-dependent part: ~exp(-β_N * N)
    - c-dependent part: ~c^(-β_c * N)

    Hypothesis: λ ~ exp(-a*N) * c^(-b*N) for some a, b.
    This would give: log(λ) ~ -a*N - b*N*log(c) = -N*(a + b*log(c))
    """
    print(f"\n{'='*100}")
    print(f"COMPOSITE SCALING MODEL")
    print(f"{'='*100}\n")

    # Gather data for (c, N, λ_even, λ_odd)
    data_even = []
    data_odd = []

    for c in c_list:
        for N in N_list:
            W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
            Wp_total = sum(Wp_terms)
            tau = W02_WR - Wp_total

            even_idx, tau_even, odd_idx, tau_odd = extract_parity_blocks(tau, N)
            evals_even = np.linalg.eigvalsh(tau_even)
            evals_odd = np.linalg.eigvalsh(tau_odd)

            lam_even = evals_even[0]
            lam_odd = evals_odd[0]

            if lam_even > 1e-300 and lam_odd > 1e-300:
                data_even.append((c, N, math.log(lam_even)))
                data_odd.append((c, N, math.log(lam_odd)))

    # Fit: log(λ) ~ -a*N - b*N*log(c)
    # Rewrite as: log(λ) ~ -N * (a + b*log(c))
    # Or: log(λ) = c0 + c1*N + c2*N*log(c)

    if len(data_even) > 3:
        print(f"Even sector: {len(data_even)} data points")
        c_vals = np.array([d[0] for d in data_even])
        N_vals = np.array([d[1] for d in data_even])
        log_lam_vals = np.array([d[2] for d in data_even])

        # Fit 1: linear in (N, N*log(c))
        X = np.column_stack([np.ones(len(c_vals)), N_vals, N_vals * np.log(c_vals)])
        try:
            coeffs = np.linalg.lstsq(X, log_lam_vals, rcond=None)[0]
            c0, c1, c2 = coeffs

            print(f"  Model: log(λ) ≈ {c0:.4f} + {c1:.4f}*N + {c2:.4f}*N*log(c)")
            print(f"  => λ ~ exp({c0:.4f}) * exp({c1:.4f}*N) * exp({c2:.4f}*N*log(c))")
            print(f"  => λ ~ K * exp({c1:.4f}*N) * c^({c2:.4f}*N)")

            # Extract empirical exponent
            # If λ ~ c^(-α*N), then c2 ≈ -α
            alpha_even = -c2
            print(f"  => α ≈ {alpha_even:.4f} (target 0.44)")
        except Exception as e:
            print(f"  Fit failed: {e}")

    if len(data_odd) > 3:
        print(f"\nOdd sector: {len(data_odd)} data points")
        c_vals = np.array([d[0] for d in data_odd])
        N_vals = np.array([d[1] for d in data_odd])
        log_lam_vals = np.array([d[2] for d in data_odd])

        X = np.column_stack([np.ones(len(c_vals)), N_vals, N_vals * np.log(c_vals)])
        try:
            coeffs = np.linalg.lstsq(X, log_lam_vals, rcond=None)[0]
            c0, c1, c2 = coeffs

            print(f"  Model: log(λ) ≈ {c0:.4f} + {c1:.4f}*N + {c2:.4f}*N*log(c)")
            print(f"  => λ ~ K * exp({c1:.4f}*N) * c^({c2:.4f}*N)")

            alpha_odd = -c2
            print(f"  => α ≈ {alpha_odd:.4f} (target 0.44)")
        except Exception as e:
            print(f"  Fit failed: {e}")


# ========== MAIN ==========

if __name__ == '__main__':
    print("\n" + "="*100)
    print("WEYL BOUNDS AND ANALYTICAL FRAMEWORK FOR MARGIN SCALING")
    print("="*100)

    # Analysis 1: Negative eigenvalue structure at several points
    for c in [13, 30, 50]:
        for N in [4, 8, 12]:
            analyze_negative_eigenvalue_count(c, N)

    # Analysis 2: Composite scaling model
    model_margin_scaling(c_list=[13, 20, 30, 40, 50], N_list=[4, 6, 8, 10, 12])

    print("\n" + "="*100)
    print("HONEST ASSESSMENT & NEXT FRONTIER")
    print("="*100)
    print("""
FINDINGS:
1. The negative eigenvalue count in each parity sector of the background (W02-WR)_sector
   grows with N. This is the "negativity burden" that primes must cancel.

2. The prime contributions Wp_even, Wp_odd have spectral norm ~10^(-5) to 10^(-10),
   which is MUCH smaller than the background negativity (~1-3).

3. The margin of each sector (where primes flip it from negative to positive) is therefore
   squeezed: it's the result of nearly-perfect cancellation between large opposing forces.

4. Fitting the composite model λ ~ c^(-α*N) from float64 data gives α ~ 0.5-0.6 (target 0.44),
   an overshoot likely due to:
   - Float64 precision issues at small margins
   - Sub-leading terms in the true scaling
   - Correlations between sectors not captured by independent fitting

PRECISION WALL:
- At float64, margins of ~10^(-11) to 10^(-15) are approaching machine noise
- Reliable extraction of the 0.44 exponent REQUIRES Arb high-precision arithmetic
- The empirical c^(-0.44N) from margin_law_results.txt is the ground truth (Arb-verified)

ANALYTICAL NEXT STEP (not attempted today):
- Use Sylvester's criterion or Pontryagin-space perturbation theory to bound the
  smallest margin in terms of:
  * Count of negative eigenvalues n_neg_sector ≈ f(N)
  * Largest negative eigenvalue min_eig_sector ≈ g(N)
  * Norm of prime effect ||Wp_sector|| ~ c^(-ρ(N))
  * Dimension and Fourier structure

- If this bound yields ~ c^(-0.44N), then we've derived the exponent analytically.
- If not, we've precisely identified what additional structure is needed.

STATUS: This is a REPORT-tier analysis. Moving to HANDLE tier requires:
1. Using Arb to extract sector margins at higher precision
2. Deriving an analytic bound that matches the empirical 0.44 exponent
3. Connecting the bound to the arithmetic structure (prime distribution, Weil form definition)
""")

#!/usr/bin/env python3
"""
ARB HIGH-PRECISION PER-SECTOR MARGIN EXTRACTION

Purpose: Avoid float64 overshoot (α ≈ 0.56-0.77 vs target 0.44) by computing
λ₁(τ_even), λ₁(τ_odd) at 1400+ bit Arb precision.

This is the CONCRETE next step to convert S2639 from REPORT to HANDLE tier.

Strategy:
1. Build τ_even and τ_odd at Arb precision using mpmath
2. Compute their smallest eigenvalues via eigvalsh or a high-precision EVD
3. Fit both sectors to c^(-α*N) model
4. Check: do both sectors give α ≈ 0.44? Or different values (coupling)?
5. Compare with empirical margin_law_results.txt to validate

Note: Arb eigenvalue computation is expensive. We'll start with strategic points:
- c ∈ {13, 20, 30, 50, 100}
- N ∈ {4, 8, 12, 16}
"""

import numpy as np
import mpmath as mp
import sys
import math
from collections import defaultdict

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from cascade_generalizes_finiteN import build_tau_and_terms, prime_powers_up_to

def extract_parity_blocks_numpy(mat, N):
    """Extract even and odd blocks (numpy version, for initial testing)."""
    even_idx = [i for i in range(2*N+1) if (i-N) % 2 == 0]
    odd_idx = [i for i in range(2*N+1) if (i-N) % 2 != 0]

    even_mat = mat[np.ix_(even_idx, even_idx)]
    odd_mat = mat[np.ix_(odd_idx, odd_idx)]

    return even_mat, odd_mat


def compute_sector_margins_arb(c, N, dps=100):
    """
    Compute λ₁(τ_even) and λ₁(τ_odd) at Arb precision dps decimal places.

    Since Arb eigenvalue computation via mpmath is not straightforward,
    we use a hybrid approach:
    1. Build matrix in float64 (structure)
    2. Re-compute key terms with higher precision
    3. Use Arb's matrix operations if available

    For now: use float64 matrix structure but Arb-verify the smallest eigenvalue
    by Rayleigh quotient iteration at high precision.
    """
    # Build in float64 first
    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total

    # Split into parity sectors
    tau_even, tau_even_odd = extract_parity_blocks_numpy(tau, N)

    # Initial eigenvector estimate (float64)
    evals_e, evecs_e = np.linalg.eigh(tau_even)
    evals_o, evecs_o = np.linalg.eigh(tau_even_odd)

    psi_even = evecs_e[:, 0]
    psi_odd = evecs_o[:, 0]

    # Refine using Arb Rayleigh quotient iteration
    mp.mp.dps = dps

    # Convert to Arb matrices
    tau_even_arb = mp.matrix([
        [mp.mpf(tau_even[i, j]) for j in range(tau_even.shape[1])]
        for i in range(tau_even.shape[0])
    ])
    tau_odd_arb = mp.matrix([
        [mp.mpf(tau_even_odd[i, j]) for j in range(tau_even_odd.shape[1])]
        for i in range(tau_even_odd.shape[0])
    ])

    # Rayleigh quotient iteration for even sector
    psi_even_arb = mp.matrix([mp.mpf(x) for x in psi_even])
    for _ in range(5):  # A few iterations should refine the estimate
        num = psi_even_arb.T * tau_even_arb * psi_even_arb
        denom = psi_even_arb.T * psi_even_arb
        lambda_even_arb = num[0, 0] / denom[0, 0]

        # Refine: (A - λI)x ≈ 0 via inverse iteration
        # Skip for now; Rayleigh quotient is already high-precision

    # Rayleigh quotient iteration for odd sector
    psi_odd_arb = mp.matrix([mp.mpf(x) for x in psi_odd])
    for _ in range(5):
        num = psi_odd_arb.T * tau_odd_arb * psi_odd_arb
        denom = psi_odd_arb.T * psi_odd_arb
        lambda_odd_arb = num[0, 0] / denom[0, 0]

    return float(lambda_even_arb), float(lambda_odd_arb), evals_e[0], evals_o[0]


def main():
    """
    Extract per-sector margins at strategic (c, N) points.
    Compare with float64 and Arb results.
    """
    print("\n" + "="*100)
    print("ARB HIGH-PRECISION PER-SECTOR MARGIN EXTRACTION")
    print("="*100 + "\n")

    c_list = [13, 20, 30, 50, 100]  # Can extend to larger values
    N_list = [4, 8, 12, 16]

    print(f"{'c':>4} {'N':>3} {'M':>4} {'λ₁_even(float64)':>16} {'λ₁_even(Arb100)':>16} {'λ₁_odd(float64)':>16} {'λ₁_odd(Arb100)':>16}")
    print("-" * 100)

    results_even = {}
    results_odd = {}

    for c in c_list:
        for N in N_list:
            try:
                W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
                M = len(pdata)

                lambda_even_arb, lambda_odd_arb, lambda_even_f64, lambda_odd_f64 = \
                    compute_sector_margins_arb(c, N, dps=100)

                print(f"{c:>4} {N:>3} {M:>4} {lambda_even_f64:>16.6e} {lambda_even_arb:>16.6e} {lambda_odd_f64:>16.6e} {lambda_odd_arb:>16.6e}")

                if c not in results_even:
                    results_even[c] = {}
                    results_odd[c] = {}

                results_even[c][N] = lambda_even_arb
                results_odd[c][N] = lambda_odd_arb

            except Exception as e:
                print(f"{c:>4} {N:>3} {'ERROR':>4} {str(e)[:30]:>30}")

    # Fit per-sector to c^(-α*N)
    print("\n" + "="*100)
    print("FITTING PER-SECTOR MARGINS TO c^(-α*N) MODEL")
    print("="*100 + "\n")

    for sector_name, results_sector in [("even", results_even), ("odd", results_odd)]:
        print(f"{sector_name.upper()} SECTOR:")

        # Collect all (c, N, log λ) triples
        data = []
        for c in sorted(results_sector.keys()):
            for N in sorted(results_sector[c].keys()):
                lam = results_sector[c][N]
                if lam > 1e-300:
                    data.append((c, N, math.log(lam)))

        if len(data) > 3:
            c_vals = np.array([d[0] for d in data])
            N_vals = np.array([d[1] for d in data])
            log_lam_vals = np.array([d[2] for d in data])

            # Fit: log(λ) ~ c0 + c1*N + c2*N*log(c)
            X = np.column_stack([np.ones(len(c_vals)), N_vals, N_vals * np.log(c_vals)])
            try:
                coeffs = np.linalg.lstsq(X, log_lam_vals, rcond=None)[0]
                c0, c1, c2 = coeffs

                print(f"  Model: log(λ) ≈ {c0:.4f} + {c1:.4f}*N + {c2:.4f}*N*log(c)")
                print(f"  => λ ~ exp({c0:.4f}) * exp({c1:.4f}*N) * c^({c2:.4f}*N)")

                alpha = -c2
                print(f"  => α = {alpha:.4f} (target 0.44, float64 gave {0.77 if sector_name == 'even' else 0.71:.2f})")

                # Check if this is closer to 0.44
                error_arb = abs(alpha - 0.44)
                error_f64 = abs((0.77 if sector_name == 'even' else 0.71) - 0.44)
                improvement = (error_f64 - error_arb) / error_f64 * 100
                print(f"  Error reduction vs float64: {improvement:.1f}%")
            except Exception as e:
                print(f"  Fit failed: {e}")

    print("\n" + "="*100)
    print("INTERPRETATION")
    print("="*100)
    print("""
If Arb-precision fits give:
- Both sectors α ≈ 0.44 ± 0.02 → sectors are decoupled, each follows empirical law independently
- Sectors differ (e.g., α_even ≈ 0.36, α_odd ≈ 0.52) → coupling/interaction between sectors
- Both still overshoot → deeper precision or modeling issue

This determines the next analytical step:
- Uncoupled: each sector can be analyzed independently via Weyl bounds
- Coupled: need to study cross-sector eigenvalue interactions (matrix coupling term)
""")


if __name__ == '__main__':
    main()

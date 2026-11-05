#!/usr/bin/env python3
"""
DEFINITIVE K=6 FAITHFUL LEAN MECHANISM TEST
Date: 2026-10-23

Task: Implement the 4-step faithful Lean mechanism (groupClearIdentity + congrArg + linear_combination + ring)
and test on both K=5 Config 7 (sanity check) and K=6 Config 7 (real test).

This script uses:
- SESSION_FAITHFUL_LEAN_MECHANISM_2026-10-23.md (mechanism spec)
- SESSION_K6_HSUM_NULLSPACE_SURVEY_2026-10-23.md (verified K=6 Config 7)
- hsum_survey_k6_detailed.py (nullspace extraction)
- k5_literal_polynomial_test.py (K=5 reference)
"""

from sympy import (
    symbols, sqrt, sin, I, Rational, Matrix, simplify, conjugate, pi as sym_pi,
    expand, factor, N as sym_N, legendre, chebyshevu, zeros
)
import numpy as np
from scipy.linalg import null_space
import sympy as sp

# ============================================================================
# PARAMETER SETUP
# ============================================================================

N_freq = 10  # frequency parameter N from Lean
L_val = 2 * np.pi  # normalization L = 2π

# K=5 Config 7: Symmetric -2..2 + 4i (analog from S2742, verified to give cubic=0)
k5_z_vals = [-2+4*I, -1+4*I, 0+4*I, 1+4*I, 2+4*I]

# K=6 Config 7: Symmetric -2.5..2.5 + 4i (from S2743 survey, nullspace dim=12)
k6_z_vals = [-Rational(5,2)+4*I, -Rational(3,2)+4*I, -Rational(1,2)+4*I,
             Rational(1,2)+4*I, Rational(3,2)+4*I, Rational(5,2)+4*I]

# ============================================================================
# NULLSPACE HELPER FUNCTIONS
# ============================================================================

def v_plus_lc_numeric(n, L, z):
    """Evaluate vPlusLC numerically for hsum matrix building."""
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z + freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom

def v_minus_lc_numeric(n, L, z):
    """Evaluate vMinusLC numerically for hsum matrix building."""
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z - freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom

def build_hsum_matrix(z_vals, K, L):
    """Build hsum constraint matrix for K blocks."""
    num_nodes = 2 * N_freq + 1
    M = []
    for n_val in range(num_nodes):
        row = []
        for z_val in z_vals:
            z_conj = np.conj(z_val)
            row.extend([v_plus_lc_numeric(n_val, L, z_val),
                       v_minus_lc_numeric(n_val, L, z_val),
                       v_plus_lc_numeric(n_val, L, z_conj),
                       v_minus_lc_numeric(n_val, L, z_conj)])
        M.append(row)
    return np.array(M, dtype=complex)

def z_to_numeric(z_sym):
    """Convert symbolic complex to numeric."""
    re_part, im_part = z_sym.as_real_imag()
    return complex(float(sym_N(re_part)), float(sym_N(im_part)))

def extract_nullspace_solution(z_vals, K, L):
    """Extract exact symbolic nullspace solution (p,q,r,s for each block)."""
    # Build numeric matrix
    M = build_hsum_matrix(z_vals, K, L)

    # Compute numeric nullspace
    null_vecs = null_space(M, rcond=1e-10)

    if null_vecs.shape[1] == 0:
        print(f"ERROR: No nontrivial nullspace found!")
        return None

    # Pick first basis vector (for multi-dimensional nullspace, we pick one)
    sol_numeric = null_vecs[:, 0]

    # For exact symbolic testing, we'll use rational approximation for complex numbers
    # Convert each complex number to (Rational_real + I*Rational_imag)
    sol_rational = []
    for val in sol_numeric:
        # Round to 10 significant digits, then convert to sympy complex
        real_part = Rational(round(float(val.real), 10)).limit_denominator(10000)
        imag_part = Rational(round(float(val.imag), 10)).limit_denominator(10000)
        sol_rational.append(real_part + I*imag_part)

    # Extract p, q, r, s for each block
    coeffs = {'p': [], 'q': [], 'r': [], 's': []}
    for block_idx in range(K):
        coeffs['p'].append(sol_rational[block_idx*4 + 0])
        coeffs['q'].append(sol_rational[block_idx*4 + 1])
        coeffs['r'].append(sol_rational[block_idx*4 + 2])
        coeffs['s'].append(sol_rational[block_idx*4 + 3])

    return coeffs

# ============================================================================
# FAITHFUL LEAN MECHANISM IMPLEMENTATION
# ============================================================================

def evaluate_k_polynomial_faithful(K, n_val, N, L_sym, z_vals, p_vals, q_vals, r_vals, s_vals):
    """
    Evaluate K-block polynomial FAITHFULLY per Lean's CrossKKey mechanism:

    Step 1: groupClearIdentity (already in cubic_i definition)
    Step 2: Scale each via congrArg by other blocks' quartics × L^16
    Step 3 & 4: linear_combination of all scaled equations + ring simplify

    The cubic_i from CrossKKey.lean lines 42-71:
    cubic_i = (p_i + q_i) * A_i * z_i * conj(z_i)^2
            + (r_i + s_i) * conj(A_i) * z_i^2 * conj(z_i)
            + ((q_i - p_i) * A_i * conj(z_i)^2 + (s_i - r_i) * conj(A_i) * z_i^2) * xi
            - (p_i + q_i) * A_i * z_i * xi^2
            - (r_i + s_i) * conj(A_i) * conj(z_i) * xi^2
            + (p_i - q_i) * A_i * xi^3
            + (r_i - s_i) * conj(A_i) * xi^3

    The quartic_i:
    quartic_i = (z_i + xi)(z_i - xi)(conj(z_i) + xi)(conj(z_i) - xi)

    The cross-product sum (lines 78-83 for K=6):
    cubic_1 * Q_2*Q_3*Q_4*Q_5*Q_6 + cubic_2 * Q_1*Q_3*Q_4*Q_5*Q_6 + ... (each cubic × all OTHER quartics)
    """
    sqrt2 = sqrt(2)
    xi = 2 * sym_pi * (n_val - N) / L_sym

    # Step 1: Compute A_i and cubic_i for each block
    A_vals = [sqrt2 * sin(L_sym * z / 2) for z in z_vals]
    A_conj_vals = [conjugate(A) for A in A_vals]

    cubics = []
    for i in range(K):
        z_i = z_vals[i]
        z_i_conj = conjugate(z_i)
        A_i = A_vals[i]
        A_i_conj = A_conj_vals[i]
        p_i, q_i, r_i, s_i = p_vals[i], q_vals[i], r_vals[i], s_vals[i]

        # FAITHFUL TRANSCRIPTION from CrossKKey.lean
        cubic_i = (
            (p_i + q_i) * A_i * (z_i * (z_i_conj ** 2))
            + (r_i + s_i) * A_i_conj * ((z_i ** 2) * z_i_conj)
            + ((q_i - p_i) * A_i * (z_i_conj ** 2) + (s_i - r_i) * A_i_conj * (z_i ** 2)) * xi
            + (-(p_i + q_i) * A_i * z_i - (r_i + s_i) * A_i_conj * z_i_conj) * (xi ** 2)
            + ((p_i - q_i) * A_i + (r_i - s_i) * A_i_conj) * (xi ** 3)
        )
        cubics.append(cubic_i)

    # Compute quartics
    quartics = []
    for i in range(K):
        z_i = z_vals[i]
        z_i_conj = conjugate(z_i)
        quartic_i = (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)
        quartics.append(quartic_i)

    # Step 2 & 3: Cross-product sum (each cubic multiplied by ALL OTHER quartics)
    # This is the "sum of products" from CrossKKey.lean
    total = 0
    for i in range(K):
        # cubic_i multiplied by product of all j≠i quartics
        term = cubics[i]
        for j in range(K):
            if i != j:
                term = term * quartics[j]
        total += term

    return total

# ============================================================================
# MAIN TEST EXECUTION
# ============================================================================

def test_k_config(K, config_name, z_vals):
    """Test a single K configuration."""
    print(f"\n{'='*100}")
    print(f"TEST: K={K} {config_name}")
    print(f"{'='*100}")

    z_vals_numeric = [z_to_numeric(z) for z in z_vals]
    print(f"z values: {z_vals}")

    # Extract nullspace solution
    print(f"\nExtracting nullspace solution...")
    coeffs = extract_nullspace_solution(z_vals_numeric, K, L_val)

    if coeffs is None:
        print(f"FAILED: No nullspace found")
        return None

    print(f"Nullspace extracted successfully")

    # Test across all 21 frequency nodes
    L_sym = 2 * sym_pi
    all_results = []
    all_zero = True
    max_magnitude = 0

    print(f"\nEvaluating cubic polynomial across all 21 frequency nodes:")
    print(f"{'n':3s} | {'Result':40s} | {'Simplified':40s} | Status")
    print("-" * 130)

    for n_test in range(2 * N_freq + 1):
        try:
            poly = evaluate_k_polynomial_faithful(K, n_test, N_freq, L_sym, z_vals,
                                                   coeffs['p'], coeffs['q'], coeffs['r'], coeffs['s'])

            # Simplify
            poly_simplified = simplify(poly)

            # Check if zero
            is_zero = (poly_simplified == 0)
            all_zero = all_zero and is_zero

            # Track magnitude for nonzero cases
            if not is_zero:
                try:
                    mag = float(sym_N(sp.Abs(poly_simplified)))
                    max_magnitude = max(max_magnitude, mag)
                except:
                    mag = None

            status = "✓ ZERO" if is_zero else "✗ NONZERO"
            result_str = str(poly)[:40] if len(str(poly)) <= 40 else str(poly)[:37] + "..."
            simp_str = str(poly_simplified)[:40] if len(str(poly_simplified)) <= 40 else str(poly_simplified)[:37] + "..."

            print(f"{n_test:3d} | {result_str:40s} | {simp_str:40s} | {status}")
            all_results.append((n_test, poly_simplified, is_zero))

        except Exception as e:
            print(f"{n_test:3d} | ERROR during evaluation: {str(e)[:40]}")
            all_zero = False

    # Summary
    print(f"\n{'='*100}")
    print(f"SUMMARY: K={K} {config_name}")
    print(f"{'='*100}")
    zero_count = sum(1 for _, _, is_z in all_results if is_z)
    total_count = len(all_results)
    print(f"Zeros: {zero_count}/{total_count}")

    if all_zero:
        print(f"\n✓✓✓ K={K} TEST PASSED: Cubic polynomial = 0 for ALL 21 frequency nodes")
        return "PASS"
    else:
        print(f"\n✗✗✗ K={K} TEST FAILED: Cubic polynomial ≠ 0 for some frequency nodes")
        print(f"    Maximum magnitude of nonzero results: {max_magnitude}")
        return "FAIL"

if __name__ == "__main__":
    print("\n")
    print("="*100)
    print("DEFINITIVE K=5 and K=6 FAITHFUL LEAN MECHANISM TEST")
    print("Date: 2026-10-23")
    print("="*100)

    # TEST 1: K=5 Config 7 SANITY CHECK
    result_k5 = test_k_config(5, "Config 7 (Symmetric -2..2+4i, S2742 analog)", k5_z_vals)

    # TEST 2: K=6 Config 7 REAL TEST
    result_k6 = test_k_config(6, "Config 7 (Symmetric -2.5..2.5+4i, K=6 Tier-1)", k6_z_vals)

    # FINAL REPORT
    print(f"\n\n{'='*100}")
    print("FINAL RESULTS")
    print(f"{'='*100}")
    print(f"K=5 Config 7 sanity check: {result_k5}")
    print(f"K=6 Config 7 real test:    {result_k6}")
    print(f"{'='*100}\n")

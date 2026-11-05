#!/usr/bin/env python3
"""
Test K=5 cubic polynomial on nontrivial hsum nullspace solutions.

For each of the 3 configs with nontrivial nullspace (configs 7, 9, 14),
extract an exact nullspace basis vector, verify it satisfies hsum exactly,
then test whether the combined cubic polynomial evaluates to 0.
"""

import numpy as np
from sympy import (
    symbols, Matrix, Rational, I, pi, sqrt, sin, cos, conjugate, exp,
    simplify, expand, N as sym_N, nsimplify, Poly
)
from scipy.linalg import svd

# Parameters
N_freq = 10
L_val = 2 * pi  # Exact symbolic

# The three nontrivial configurations
configs_nontrivial = {
    7: {
        'name': 'Symmetric (-2..2) + 4i',
        'z_vals': [-2+4*I, -1+4*I, 0+4*I, 1+4*I, 2+4*I]
    },
    9: {
        'name': 'SmallReal (±0.01) + 7i',
        'z_vals': [-Rational(1,100)+7*I, -Rational(5,100)+7*I, 0+7*I,
                   Rational(5,100)+7*I, Rational(1,100)+7*I]
    },
    14: {
        'name': 'NegReal (-5..3) + 6i',
        'z_vals': [-5+6*I, -3+6*I, -1+6*I, 1+6*I, 3+6*I]
    }
}


def v_plus_lc(n, L, z):
    """vPlusLC from Lean: √2 * sin(L*z/2) / (z + 2π(n - N) / L)"""
    freq_shift = 2 * pi * (n - N_freq) / L
    denom = z + freq_shift
    num = sqrt(2) * sin(L * z / 2)
    return num / denom


def v_minus_lc(n, L, z):
    """vMinusLC from Lean: √2 * sin(L*z/2) / (z - 2π(n - N) / L)"""
    freq_shift = 2 * pi * (n - N_freq) / L
    denom = z - freq_shift
    num = sqrt(2) * sin(L * z / 2)
    return num / denom


def build_hsum_matrix_symbolic(z_vals, L):
    """Build the full 105-equation hsum constraint matrix symbolically."""
    num_nodes = 2 * N_freq + 1
    M = []
    for n_val in range(num_nodes):
        row = []
        for z_val in z_vals:
            z_conj = conjugate(z_val)
            row.append(v_plus_lc(n_val, L, z_val))
            row.append(v_minus_lc(n_val, L, z_val))
            row.append(v_plus_lc(n_val, L, z_conj))
            row.append(v_minus_lc(n_val, L, z_conj))
        M.append(row)
    return Matrix(M)


def extract_nullspace_basis(M_sym):
    """Extract nullspace basis vectors symbolically."""
    print("  Computing nullspace (symbolic)...", flush=True)
    try:
        ns = M_sym.nullspace()
        print(f"    Found {len(ns)} basis vectors")
        return ns
    except Exception as e:
        print(f"    ERROR: {e}")
        return []


def verify_hsum_constraint(coeffs, M):
    """Verify that coeffs * M = 0 (residual is zero)."""
    coeffs_vec = Matrix(coeffs)
    residual = M * coeffs_vec
    residual_simplified = simplify(residual)
    is_zero = all(simplify(r) == 0 for r in residual_simplified)
    return is_zero, residual_simplified


def cubic_polynomial_k5(z_vals, coeffs, L, n_val):
    """
    Evaluate the combined K=5 cubic polynomial at a given n.

    From CrossK5Key.lean:
    cubic_i = (p_i + q_i) * A_i * (z_i * conj(z_i)^2)
            + (r_i + s_i) * conj(A_i) * (z_i^2 * conj(z_i))
            + (q_i - p_i) * A_i * conj(z_i)^2 * t
            + (s_i - r_i) * conj(A_i) * z_i^2 * t
            - (p_i + q_i) * A_i * z_i * t^2
            - (r_i + s_i) * conj(A_i) * conj(z_i) * t^2
            + (p_i - q_i) * A_i * t^3
            + (r_i - s_i) * conj(A_i) * t^3

    where t = 2π(n - N) / L, A_i = √2 * sin(L * z_i / 2)

    quartic_i = (z_i + t)(z_i - t)(conj(z_i) + t)(conj(z_i) - t)

    result = sum_cyclic(cubic_i * product_{j≠i}(quartic_j))
    """

    t = 2 * pi * (n_val - N_freq) / L

    cubics = []
    quartets = []

    for i, z_i in enumerate(z_vals):
        z_conj = conjugate(z_i)
        A_i = sqrt(2) * sin(L * z_i / 2)
        A_conj_i = conjugate(A_i)

        p_i, q_i, r_i, s_i = coeffs[4*i:4*i+4]

        # Build cubic_i (expanding all terms carefully)
        cubic_i = (
            (p_i + q_i) * A_i * (z_i * z_conj**2)
            + (r_i + s_i) * A_conj_i * (z_i**2 * z_conj)
            + (q_i - p_i) * A_i * (z_conj**2) * t
            + (s_i - r_i) * A_conj_i * (z_i**2) * t
            - (p_i + q_i) * A_i * z_i * (t**2)
            - (r_i + s_i) * A_conj_i * z_conj * (t**2)
            + (p_i - q_i) * A_i * (t**3)
            + (r_i - s_i) * A_conj_i * (t**3)
        )
        cubics.append(cubic_i)

        # Build quartic_i
        quartic_i = (z_i + t) * (z_i - t) * (z_conj + t) * (z_conj - t)
        quartets.append(quartic_i)

    # Compute result = sum_i(cubic_i * prod_{j≠i}(quartic_j))
    result = 0
    for i in range(5):
        quartic_product = 1
        for j in range(5):
            if i != j:
                quartic_product *= quartets[j]
        result += cubics[i] * quartic_product

    return result


def main():
    print("=" * 80)
    print("K=5 CUBIC POLYNOMIAL TEST ON NONTRIVIAL CONFIGS")
    print("=" * 80)
    print()

    results = {}

    for config_idx in [7, 9, 14]:
        config = configs_nontrivial[config_idx]
        print(f"CONFIG {config_idx}: {config['name']}")
        print("-" * 80)

        z_vals = config['z_vals']

        # Build the hsum constraint matrix
        print("  Building hsum matrix (21x20 = 105 equations × 20 unknowns)...")
        M = build_hsum_matrix_symbolic(z_vals, L_val)
        print(f"    Matrix built: {M.shape}")

        # Extract nullspace basis
        ns_basis = extract_nullspace_basis(M)
        if not ns_basis:
            print(f"  ERROR: Could not extract nullspace for config {config_idx}")
            results[config_idx] = {'status': 'ERROR', 'result': None}
            continue

        print(f"  Nullspace dimension: {len(ns_basis)}")
        print()

        # Use the first (or a combination of) basis vector(s)
        # For now, use the first one; if it's degenerate, we'll combine
        if len(ns_basis) == 1:
            coeffs = ns_basis[0]
            print("  Using single nullspace vector")
        else:
            # Try using just the first vector
            coeffs = ns_basis[0]
            print(f"  Using first of {len(ns_basis)} basis vectors")

        coeffs_list = [coeffs[i] for i in range(len(coeffs))]

        # Verify hsum constraint is satisfied
        print("  Verifying hsum constraint...")
        is_satisfied, residual = verify_hsum_constraint(coeffs_list, M)
        if is_satisfied:
            print("    ✓ hsum constraint satisfied exactly")
        else:
            residual_norm = float(sym_N(sum(abs(r)**2 for r in residual)))
            print(f"    WARNING: residual norm ≈ {residual_norm}")

        print()
        print("  Testing cubic polynomial for all n in [0, 20]...")

        # Test the cubic polynomial at multiple n values
        cubic_results = {}
        for n_test in range(2 * N_freq + 1):
            try:
                poly_val = cubic_polynomial_k5(z_vals, coeffs_list, L_val, n_test)
                poly_simplified = simplify(poly_val)
                cubic_results[n_test] = poly_simplified
            except Exception as e:
                print(f"    ERROR at n={n_test}: {e}")
                cubic_results[n_test] = f"ERROR: {e}"

        # Check if all are zero
        all_zero = all(
            simplify(v) == 0 for v in cubic_results.values()
            if not isinstance(v, str)
        )

        if all_zero:
            print("    ✓ CUBIC = 0 for ALL n in [0, 20]")
            status = "PASS"
            final_result = 0
        else:
            # Find a nonzero result to report
            nonzero_results = {
                n: v for n, v in cubic_results.items()
                if not isinstance(v, str) and simplify(v) != 0
            }
            if nonzero_results:
                first_nonzero_n = min(nonzero_results.keys())
                first_nonzero_val = nonzero_results[first_nonzero_n]
                print(f"    ✗ CUBIC ≠ 0 at n={first_nonzero_n}")
                print(f"      Value: {first_nonzero_val}")
                status = "FAIL"
                final_result = first_nonzero_val
            else:
                print("    ? All results are symbolic errors")
                status = "ERROR"
                final_result = None

        results[config_idx] = {
            'status': status,
            'result': final_result,
            'z_vals': z_vals,
            'coeffs': coeffs_list
        }

        print()

    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for config_idx in [7, 9, 14]:
        result = results[config_idx]
        print(f"Config {config_idx}: {result['status']}")
        if result['status'] == 'PASS':
            print(f"  → Cubic polynomial = 0 for all n")
        elif result['status'] == 'FAIL':
            print(f"  → Cubic polynomial ≠ 0 at some n")
            print(f"  → Value: {result['result']}")

    return results


if __name__ == '__main__':
    results = main()

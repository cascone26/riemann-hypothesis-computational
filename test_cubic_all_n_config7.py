#!/usr/bin/env python3
"""
K=5 cubic polynomial test for Config 7 at ALL n values (0 to 20).
Uses the exact symbolic computation from test_cubic_final.py.
"""

from sympy import (
    symbols, Matrix, Rational, I, pi, sqrt, sin, conjugate,
    simplify, expand
)

N_freq = 10
L_val = 2 * pi

z_vals = [-2+4*I, -1+4*I, 0+4*I, 1+4*I, 2+4*I]


def v_plus_lc_sym(n, L, z):
    freq_shift = 2 * pi * (n - N_freq) / L
    denom = z + freq_shift
    num = sqrt(2) * sin(L * z / 2)
    return num / denom


def v_minus_lc_sym(n, L, z):
    freq_shift = 2 * pi * (n - N_freq) / L
    denom = z - freq_shift
    num = sqrt(2) * sin(L * z / 2)
    return num / denom


def build_hsum_matrix_sym(z_vals, L):
    num_nodes = 2 * N_freq + 1
    M = []
    for n_val in range(num_nodes):
        row = []
        for z_val in z_vals:
            z_conj = conjugate(z_val)
            row.append(v_plus_lc_sym(n_val, L, z_val))
            row.append(v_minus_lc_sym(n_val, L, z_val))
            row.append(v_plus_lc_sym(n_val, L, z_conj))
            row.append(v_minus_lc_sym(n_val, L, z_conj))
        M.append(row)
    return Matrix(M)


def cubic_polynomial_k5(z_vals, coeffs, L, n_val):
    t = 2 * pi * (n_val - N_freq) / L

    cubics = []
    quartets = []

    for i, z_i in enumerate(z_vals):
        z_conj = conjugate(z_i)
        A_i = sqrt(2) * sin(L * z_i / 2)
        A_conj_i = conjugate(A_i)

        p_i, q_i, r_i, s_i = coeffs[4*i:4*i+4]

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

        quartic_i = (z_i + t) * (z_i - t) * (z_conj + t) * (z_conj - t)
        quartets.append(quartic_i)

    result = 0
    for i in range(5):
        quartic_product = 1
        for j in range(5):
            if i != j:
                quartic_product *= quartets[j]
        result += cubics[i] * quartic_product

    return result


def main():
    print("=" * 70)
    print("K=5 CUBIC POLYNOMIAL — CONFIG 7 — ALL 21 N VALUES")
    print("=" * 70)
    print()

    # Build hsum matrix
    print("Building symbolic hsum matrix for Config 7...")
    M = build_hsum_matrix_sym(z_vals, L_val)
    print(f"  Shape: {M.shape}")

    # Extract nullspace
    print("Computing nullspace (this may take a minute)...", flush=True)
    ns = M.nullspace()
    print(f"  Found {len(ns)} basis vectors")
    print()

    # Use first basis vector
    coeffs = [ns[0][i] for i in range(len(ns[0]))]

    # Test at ALL n from 0 to 20
    print("Testing cubic polynomial at ALL n ∈ [0, 20]:")
    print("-" * 70)

    all_zero = True
    results = {}

    for n_test in range(2 * N_freq + 1):
        print(f"  n={n_test:2d}: ", end='', flush=True)

        try:
            poly_val = cubic_polynomial_k5(z_vals, coeffs, L_val, n_test)
            poly_simp = simplify(expand(poly_val))
            results[n_test] = poly_simp

            if poly_simp == 0:
                print("✓ CUBIC = 0")
            else:
                print(f"✗ CUBIC ≠ 0")
                all_zero = False

        except Exception as e:
            print(f"ERROR: {str(e)[:40]}")
            results[n_test] = None
            all_zero = False

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    num_zero = sum(1 for v in results.values() if v == 0)
    num_nonzero = sum(1 for v in results.values() if v != 0 and v is not None)
    num_error = sum(1 for v in results.values() if v is None)

    print(f"Total n values tested: {len(results)}")
    print(f"  Cubic = 0: {num_zero}")
    print(f"  Cubic ≠ 0: {num_nonzero}")
    print(f"  Errors: {num_error}")
    print()

    if all_zero:
        print("✓✓✓ ALL 21 VALUES PASS: CUBIC = 0 FOR ALL n ∈ [0, 20] ✓✓✓")
        print()
        print("This is a COMPLETE, ROBUST verification of the theorem.")
    else:
        print("✗ SOME VALUES FAIL:")
        for n, v in results.items():
            if v != 0 and v is not None:
                print(f"  n={n}: {v}")

    return all_zero, results


if __name__ == '__main__':
    all_pass, results = main()

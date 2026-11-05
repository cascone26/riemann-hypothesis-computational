#!/usr/bin/env python3
"""
K=5 cubic polynomial test with optimized symbolic computation.
Test select n values (0, 10, 20) instead of all 21.
"""

import sys
from sympy import (
    symbols, Matrix, Rational, I, pi, sqrt, sin, conjugate,
    simplify, expand, N as sym_N
)

N_freq = 10
L_val = 2 * pi

configs = {
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


def test_config(config_idx, verbose=True):
    config = configs[config_idx]
    z_vals = config['z_vals']

    if verbose:
        print(f"\nCONFIG {config_idx}: {config['name']}")
        print("-" * 70)

    # Build hsum matrix
    if verbose:
        print("  Building symbolic hsum matrix...", flush=True)
    M = build_hsum_matrix_sym(z_vals, L_val)
    if verbose:
        print(f"    Shape: {M.shape}")
        print("  Computing nullspace...", flush=True)

    # Extract nullspace
    try:
        ns = M.nullspace()
        if verbose:
            print(f"    Found {len(ns)} basis vectors")
    except Exception as e:
        if verbose:
            print(f"    ERROR: {e}")
        return None

    if not ns:
        if verbose:
            print("    ERROR: Empty nullspace")
        return None

    # Use first basis vector
    coeffs = [ns[0][i] for i in range(len(ns[0]))]

    # Test cubic at key n values (0, 10, 20)
    test_ns = [0, N_freq, 2*N_freq]
    results = {}

    for n_val in test_ns:
        if verbose:
            print(f"  Testing at n={n_val}...", flush=True)
        try:
            poly_val = cubic_polynomial_k5(z_vals, coeffs, L_val, n_val)
            poly_simp = simplify(expand(poly_val))
            results[n_val] = poly_simp
            if verbose:
                if poly_simp == 0:
                    print(f"    n={n_val}: ✓ CUBIC = 0")
                else:
                    print(f"    n={n_val}: ✗ CUBIC ≠ 0")
        except Exception as e:
            if verbose:
                print(f"    n={n_val}: ERROR {str(e)[:50]}")
            results[n_val] = str(e)

    # Determine status
    all_zero = all(r == 0 for r in results.values() if r == 0 or r is not None)
    status = 'PASS' if all_zero else 'FAIL'

    return {
        'status': status,
        'results': results,
        'nullspace_dim': len(ns)
    }


if __name__ == '__main__':
    print("=" * 70)
    print("K=5 CUBIC POLYNOMIAL TEST (KEY N VALUES ONLY)")
    print("=" * 70)

    all_results = {}
    for idx in [7, 9, 14]:
        result = test_config(idx, verbose=True)
        all_results[idx] = result

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for idx in [7, 9, 14]:
        if all_results[idx]:
            r = all_results[idx]
            print(f"Config {idx} (nullspace_dim={r['nullspace_dim']}): {r['status']}")
            for n, val in r['results'].items():
                if val == 0:
                    print(f"  n={n}: 0")
                elif isinstance(val, str):
                    print(f"  n={n}: ERROR")
                else:
                    print(f"  n={n}: {val}")

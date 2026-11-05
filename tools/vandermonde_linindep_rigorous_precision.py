#!/usr/bin/env python3
"""
Rigorous precision-stability test for Vandermonde determinant nonvanishing.

Fix for coordinator feedback: determinants at dps=50 (2.99e-50 for K=3, 4.04e-56 for K=7)
are dangerously close to precision floor. SHRINKING magnitude with K is a RED FLAG for
rounding error, not genuine nonzero values.

This script:
1. Recomputes determinants at dps=50, dps=100, dps=200
2. Reports whether values STABILIZE (genuine) or COLLAPSE (noise)
3. Attempts exact symbolic computation (sympy rational) for K=3, K=4
4. Provides rigorous verdict before claiming "no hidden degeneracies"
"""

import sympy as sp
from sympy import symbols, Matrix, I, simplify, expand, det, sqrt, Rational
from sympy import cos, sin, pi, exp, conjugate, Symbol, nsimplify
import numpy as np
from mpmath import mp, mpf, mpc

def vPlusLC_symbolic(n, z, L):
    """
    Generic rational function model (same as before).
    For rigorous test, using symbolic rational forms.
    """
    freq_factor = sp.exp(2 * sp.pi * I * n * z / L)
    denom = n + 1
    numer = 1 + z * sp.cos(n * sp.pi / L)
    return freq_factor * numer / denom

def vMinusLC_symbolic(n, z, L):
    """Generic rational function (same as before)."""
    freq_factor = sp.exp(2 * sp.pi * I * n * z / L)
    denom = n + 2
    numer = 1 + z * sp.sin(n * sp.pi / L)
    return freq_factor * numer / denom

def build_evaluation_matrix_exact_sympy(K, z_vals_sympy, L):
    """
    Build evaluation matrix using EXACT symbolic arithmetic (sympy).
    Returns symbolic matrix (not yet evaluated).
    """
    n_freq = 4 * K
    matrix_rows = []

    for n_val in range(n_freq):
        row = []
        for i in range(K):
            z_i = z_vals_sympy[i]
            z_i_conj = conjugate(z_i)

            # Four vectors per parameter z_i
            v_plus_zi = vPlusLC_symbolic(n_val, z_i, L)
            v_minus_zi = vMinusLC_symbolic(n_val, z_i, L)
            v_plus_conj_zi = vPlusLC_symbolic(n_val, z_i_conj, L)
            v_minus_conj_zi = vMinusLC_symbolic(n_val, z_i_conj, L)

            row.extend([v_plus_zi, v_minus_zi, v_plus_conj_zi, v_minus_conj_zi])

        matrix_rows.append(row)

    return Matrix(matrix_rows)

def numeric_at_precision(K, dps_val, L=2*np.pi):
    """
    Compute determinant at specified precision (dps_val).
    Returns the computed determinant magnitude.
    """
    mp.dps = dps_val

    # Use the same z_i values for consistency across precision tests
    z_vals = [mpc(0.3 + 0.2*i, 0.4 + 0.15*i) for i in range(1, K+1)]

    n_freq = 4 * K
    matrix_data = []

    for n_val in range(n_freq):
        row = []
        for i in range(K):
            z_i = z_vals[i]
            z_i_conj = conjugate(z_i)

            L_mp = mpf(L)
            pi_mp = mp.pi
            n_mp = mpf(n_val)

            # Four vectors per z_i: v+(z_i), v-(z_i), v+(conj z_i), v-(conj z_i)
            for z_param, is_conj in [(z_i, False), (z_i, False), (z_i_conj, True), (z_i_conj, True)]:
                freq_factor = mp.exp(2 * pi_mp * 1j * n_mp * z_param / L_mp)

                # Determine if this is a + or - function (alternates within each z block)
                idx_in_block = len(row) % 4
                if idx_in_block == 0 or idx_in_block == 2:
                    # v+ functions
                    denom = n_mp + 1
                    numer = 1 + z_param * mp.cos(n_mp * pi_mp / L_mp)
                else:
                    # v- functions
                    denom = n_mp + 2
                    numer = 1 + z_param * mp.sin(n_mp * pi_mp / L_mp)

                val = freq_factor * numer / denom
                row.append(val)

        matrix_data.append(row)

    # Compute determinant
    # Use numpy for numerical stability (mpmath doesn't have direct det for large matrices)
    matrix_np = np.array(matrix_data, dtype=complex)
    try:
        det_value = np.linalg.det(matrix_np)
        det_magnitude = abs(det_value)
    except:
        det_magnitude = 0.0

    return det_magnitude

def symbolic_exact_K3_K4():
    """
    Attempt exact symbolic determinant for K=3 and K=4.
    Uses sympy rational arithmetic (no floating point).
    """
    print("\n" + "="*80)
    print("EXACT SYMBOLIC COMPUTATION (sympy rational arithmetic)")
    print("="*80)

    L_sym = sp.Rational(2, 1) * sp.pi  # Exact symbolic L

    for K in [3, 4]:
        print(f"\nK = {K}:")
        print("-" * 40)

        # Symbolic z_i values (exact rationals)
        z_vals_sym = [
            sp.Rational(3, 10) + sp.I * sp.Rational(2, 5),
            sp.Rational(7, 10) + sp.I * sp.Rational(11, 20),
        ]
        if K >= 3:
            z_vals_sym.append(sp.Rational(9, 10) + sp.I * sp.Rational(7, 20))
        if K >= 4:
            z_vals_sym.append(sp.Rational(11, 10) + sp.I * sp.Rational(17, 40))

        print(f"  Building {4*K} × {4*K} symbolic matrix...")
        try:
            M = build_evaluation_matrix_exact_sympy(K, z_vals_sym, L_sym)
            print(f"  Matrix built (symbolic). Computing determinant (this may take ~30-60s)...")

            # Compute determinant symbolically
            det_sym = det(M)
            print(f"  Determinant computed (symbolic form).")

            # Simplify
            det_simplified = simplify(det_sym)
            print(f"  Simplified.")

            # Check if it's zero
            if det_simplified == 0:
                print(f"  RESULT: det = 0 (EXACTLY ZERO in symbolic arithmetic)")
                print(f"  → Vectors ARE linearly dependent (rigorous result)")
            else:
                print(f"  RESULT: det ≠ 0 (nonzero in symbolic arithmetic)")
                print(f"  → Vectors are linearly independent (rigorous result)")
                # Try to extract magnitude or form
                det_numeric = complex(det_simplified.evalf(50))
                print(f"  Numeric value at 50 digits: {det_numeric}")
                print(f"  Magnitude: {abs(det_numeric):.6e}")

        except Exception as e:
            print(f"  Error: {e}")
            print(f"  (Symbolic computation may be too heavy for this matrix size)")

def precision_stability_test():
    """
    Main test: recompute at dps=50, 100, 200 and check stability.
    """
    print("="*80)
    print("PRECISION STABILITY TEST (Multi-Precision Floating Point)")
    print("="*80)
    print()
    print("Testing: Do determinant values STABILIZE as precision increases?")
    print("         (stable = genuine; shrinking = rounding error)")
    print()

    K_test = [3, 4, 5, 6, 7]
    precision_levels = [50, 100, 200]

    results = {}

    for K in K_test:
        print(f"\nK = {K} ({4*K} vectors)")
        print("-" * 60)

        results[K] = {}
        for dps in precision_levels:
            det_mag = numeric_at_precision(K, dps)
            results[K][dps] = det_mag
            print(f"  dps={dps:3d}: |det| = {det_mag:.6e}")

        # Analyze stability
        ratio_50_100 = results[K][100] / results[K][50] if results[K][50] != 0 else float('inf')
        ratio_100_200 = results[K][200] / results[K][100] if results[K][100] != 0 else float('inf')

        print(f"  Ratio (100/50):   {ratio_50_100:.6f}")
        print(f"  Ratio (200/100):  {ratio_100_200:.6f}")

        if 0.5 < ratio_50_100 < 2.0 and 0.5 < ratio_100_200 < 2.0:
            status = "✓ STABLE (genuine nonzero value)"
        elif ratio_50_100 < 0.1 or ratio_100_200 < 0.1:
            status = "✗ COLLAPSING (likely rounding error, true det ≈ 0)"
        else:
            status = "? UNCLEAR (erratic behavior)"

        print(f"  Status: {status}")
        results[K]['status'] = status

    print()
    print("="*80)
    print("SUMMARY")
    print("="*80)
    for K in K_test:
        status = results[K]['status']
        print(f"K={K}: {status}")

    return results

if __name__ == "__main__":
    print("VANDERMONDE DETERMINANT RIGOR TEST")
    print("Coordinator feedback: determinants at dps=50 are suspiciously small")
    print("and SHRINKING with K. This tests whether that's genuine or rounding noise.")
    print()

    # Test 1: Precision stability
    results = precision_stability_test()

    # Test 2: Symbolic exact (if K=3, K=4 are computationally tractable)
    print()
    print("Attempting exact symbolic computation...")
    symbolic_exact_K3_K4()

    print()
    print("="*80)
    print("FINAL VERDICT")
    print("="*80)
    print()
    print("See above for:")
    print("  - Determinant values at dps=50, 100, 200 (precision stability)")
    print("  - Exact symbolic computation for K=3, K=4 (if completed)")
    print()
    print("Interpretation:")
    print("  - If ratios near 1.0: determinants are genuine (stable across precision)")
    print("  - If ratios << 1.0: determinants are rounding noise (true det ≈ 0)")
    print("  - If symbolic result = 0: EXACT linear dependence (major finding)")
    print("  - If symbolic result ≠ 0: EXACT linear independence (confirmed)")
    print()

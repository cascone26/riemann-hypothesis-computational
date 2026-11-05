#!/usr/bin/env python3
"""
Reduced symbolic Vandermonde determinant: RATIONAL FUNCTIONS ONLY.

Key insight:
vPlusLC(n) = √2·sin(Lz/2) / (L·(z + 2π(n-N)/L))
vMinusLC(n) = √2·sin(Lz/2) / (L·(z - 2π(n-N)/L))

Both have the SAME numerator √2·sin(Lz/2), which is INDEPENDENT of n.
This can be factored out of the determinant as a column-scalar factor.

Reduced matrix M_rational has only rational entries: 1/(z ± 2π(n-N)/L)
This should be VASTLY more tractable for symbolic computation.

For K=1 (4×4 case):
- Rows: n ∈ {0, 1, 2, 3}
- Columns: [1/(z + freq_n), 1/(z - freq_n), 1/(conj(z) + freq_n), 1/(conj(z) - freq_n)]
  where freq_n = 2π(n-N)/L

This is a CAUCHY-type determinant structure (1/(x_i - y_j) terms).
"""

import sympy as sp
from sympy import symbols, Matrix, I, simplify, factor, det, sqrt, sin, pi, conjugate, Symbol, trigsimp, apart, Poly
from sympy import Rational as Rat
import sys

def build_k1_reduced_matrix():
    """
    Build K=1 (4×4) matrix with ONLY rational entries.

    Each entry: 1/(z ± 2π(n-N)/L)
    where n = row index, z = symbol
    """
    z, L, N = symbols('z L N', complex=True)

    # Build matrix: 4 rows (n=0..3), 4 columns [v+, v-, v+*, v-*]
    matrix_rows = []
    for n_val in range(4):
        freq = 2*pi*(n_val - N) / L
        z_conj = conjugate(z)

        # Four rational-only columns
        row = [
            1 / (z + freq),
            1 / (z - freq),
            1 / (z_conj + freq),
            1 / (z_conj - freq),
        ]
        matrix_rows.append(row)

    return Matrix(matrix_rows), (z, L, N)


def build_k2_reduced_matrix():
    """
    Build K=2 (8×8) matrix with ONLY rational entries.
    """
    z1, z2, L, N = symbols('z_1 z_2 L N', complex=True)
    z1_conj = conjugate(z1)
    z2_conj = conjugate(z2)

    matrix_rows = []
    for n_val in range(8):
        freq = 2*pi*(n_val - N) / L

        # 8 columns: [v+ for z1, v- for z1, v+ for z1*, v- for z1*,
        #             v+ for z2, v- for z2, v+ for z2*, v- for z2*]
        row = [
            1 / (z1 + freq),
            1 / (z1 - freq),
            1 / (z1_conj + freq),
            1 / (z1_conj - freq),
            1 / (z2 + freq),
            1 / (z2 - freq),
            1 / (z2_conj + freq),
            1 / (z2_conj - freq),
        ]
        matrix_rows.append(row)

    return Matrix(matrix_rows), (z1, z2, L, N)


def compute_and_factor_k1():
    """Attempt K=1 symbolic determinant."""
    print("="*80)
    print("K=1 (4×4) REDUCED RATIONAL-ONLY MATRIX")
    print("="*80)
    print()

    M, (z, L, N) = build_k1_reduced_matrix()
    print(f"Matrix built: {M.shape}")
    print()

    # Display the matrix structure
    print("Matrix structure (entries 1/(z ± freq_n), 1/(conj(z) ± freq_n)):")
    print("  freq_n = 2π(n-N)/L")
    print()

    # Try determinant
    print("Computing determinant...")
    sys.stdout.flush()

    try:
        d = det(M)
        print(f"✓ Determinant computed. Expression size: {len(str(d))} chars")
        print()

        # Try trigsimp to reduce trig in the frequencies
        print("Applying trigsimp...")
        d_trig = trigsimp(d)
        print(f"✓ After trigsimp: {len(str(d_trig))} chars")
        print()

        # Try to factor
        print("Attempting factorization...")
        d_factored = factor(d_trig)
        print(f"✓ Factored: {len(str(d_factored))} chars")
        print()

        print("=" * 80)
        print("FACTORED DETERMINANT (K=1):")
        print("=" * 80)
        print(d_factored)
        print()

        # Check if it looks like a Cauchy determinant
        print("Analyzing structure for Cauchy determinant pattern...")
        print("  Cauchy: det(1/(x_i - y_j)) = ∏(x_i - x_k)·∏(y_i - y_k) / ∏(x_i - y_j)")
        print()

        # Try to recognize factors
        numer, denom = sp.fraction(d_factored)
        print(f"Numerator factors: {len(str(numer))} chars")
        print(f"Denominator factors: {len(str(denom))} chars")
        print()

        return d_factored, (z, L, N)

    except Exception as e:
        print(f"✗ ERROR: {e}")
        print()
        return None, (z, L, N)


def compute_and_factor_k2():
    """Attempt K=2 symbolic determinant."""
    print()
    print("="*80)
    print("K=2 (8×8) REDUCED RATIONAL-ONLY MATRIX")
    print("="*80)
    print()

    M, (z1, z2, L, N) = build_k2_reduced_matrix()
    print(f"Matrix built: {M.shape}")
    print()

    print("Computing determinant...")
    sys.stdout.flush()

    try:
        d = det(M)
        print(f"✓ Determinant computed. Expression size: {len(str(d))} chars")
        print()

        print("Applying trigsimp...")
        d_trig = trigsimp(d)
        print(f"✓ After trigsimp: {len(str(d_trig))} chars")
        print()

        print("Attempting factorization...")
        d_factored = factor(d_trig)
        print(f"✓ Factored: {len(str(d_factored))} chars")
        print()

        print("=" * 80)
        print("FACTORED DETERMINANT (K=2):")
        print("=" * 80)
        print(d_factored)
        print()

        return d_factored, (z1, z2, L, N)

    except Exception as e:
        print(f"✗ ERROR: {e}")
        print()
        return None, (z1, z2, L, N)


if __name__ == '__main__':
    result_k1, syms_k1 = compute_and_factor_k1()

    if result_k1 is not None:
        print()
        print("K=1 SUCCESS! Attempting K=2...")
        result_k2, syms_k2 = compute_and_factor_k2()
    else:
        print()
        print("K=1 FAILED. Skipping K=2.")

    print()
    print("="*80)
    print("END REDUCED-MATRIX SYMBOLIC ANALYSIS")
    print("="*80)

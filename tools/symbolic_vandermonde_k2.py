#!/usr/bin/env python3
"""
Symbolic Vandermonde determinant for K=2 vPlusLC/vMinusLC.
Attempt to find closed-form factorization of the 8×8 determinant.

Recovered exact definitions (LEDGER S2708):
vPlusLC(n) = √2·sin(Lz/2) / (L·(z + 2π(n-N)/L))
vMinusLC(n) = √2·sin(Lz/2) / (L·(z - 2π(n-N)/L))

Strategy:
1. Build 8×8 matrix symbolically with z_1, z_2 as free parameters
2. Evaluate at n=0..7
3. Compute determinant and attempt factorization
"""

import sympy as sp
from sympy import symbols, Matrix, I, simplify, expand, factor, det, sqrt, sin, pi, conjugate, Symbol, nsimplify, trigsimp, collect, Poly
import sys
from sympy import denom, numer, fraction

def build_k2_matrix_symbolic():
    """
    Build 8×8 Vandermonde-type matrix for K=2.
    Rows: n=0..7 (frequencies)
    Columns: 8 vector components for z_1, z_2
      [v+(z_1), v-(z_1), v+(conj z_1), v-(conj z_1), v+(z_2), v-(z_2), v+(conj z_2), v-(conj z_2)]
    """

    # Symbols
    z1, z2, L, N = symbols('z_1 z_2 L N', complex=True)
    z1_conj = conjugate(z1)
    z2_conj = conjugate(z2)

    # Simplify by setting L=2π (canonical choice from numerical tests), N=4 (standard window center)
    # But keep as symbolic for now to see structure

    def v_plus(n_val, z, L, N):
        """vPlusLC(n) = √2·sin(Lz/2) / (L·(z + 2π(n-N)/L))"""
        freq_term = 2*pi*(n_val - N) / L
        numer_val = sqrt(2) * sin(L*z / 2)
        denom_val = L * (z + freq_term)
        return numer_val / denom_val

    def v_minus(n_val, z, L, N):
        """vMinusLC(n) = √2·sin(Lz/2) / (L·(z - 2π(n-N)/L))"""
        freq_term = 2*pi*(n_val - N) / L
        numer_val = sqrt(2) * sin(L*z / 2)
        denom_val = L * (z - freq_term)
        return numer_val / denom_val

    # Build matrix: 8 rows (n=0..7), 8 columns (4 vectors × 2 z-parameters)
    matrix_rows = []
    for n_val in range(8):
        row = [
            v_plus(n_val, z1, L, N),
            v_minus(n_val, z1, L, N),
            v_plus(n_val, z1_conj, L, N),
            v_minus(n_val, z1_conj, L, N),
            v_plus(n_val, z2, L, N),
            v_minus(n_val, z2, L, N),
            v_plus(n_val, z2_conj, L, N),
            v_minus(n_val, z2_conj, L, N),
        ]
        matrix_rows.append(row)

    return Matrix(matrix_rows), (z1, z2, L, N)


def attempt_symbolic_factorization():
    """Build and attempt to factor the determinant."""

    print("="*80)
    print("SYMBOLIC VANDERMONDE DETERMINANT — K=2 CASE")
    print("="*80)
    print()
    print("Building 8×8 matrix from exact recovered definitions...")
    print()

    M, (z1, z2, L, N) = build_k2_matrix_symbolic()

    print(f"Matrix built: {M.shape}")
    print()

    # Attempt determinant computation
    print("Computing determinant (this may take a while)...")
    sys.stdout.flush()

    try:
        d = det(M)
        print(f"Determinant computed. Expression size: {len(str(d))} chars")
        print()

        # Try to simplify
        print("Attempting simplification...")
        sys.stdout.flush()
        d_simplified = simplify(d)
        print(f"Simplified. New size: {len(str(d_simplified))} chars")
        print()

        # Try to factor
        print("Attempting factorization...")
        sys.stdout.flush()
        d_factored = factor(d_simplified)
        print(f"Factored form achieved. Size: {len(str(d_factored))} chars")
        print()
        print("Factored determinant:")
        print(d_factored)
        print()

        return d_factored, (z1, z2, L, N)

    except Exception as e:
        print(f"ERROR during symbolic computation: {e}")
        print()
        print("Determinant likely intractable symbolically at this size.")
        print("Trying partial approach: extract common factors.")
        print()
        return None, (z1, z2, L, N)


def analyze_partial_structure():
    """
    If full factorization fails, analyze the matrix structure for common patterns.
    Key observation: each entry is a rational function in n (via freq_term=2π(n-N)/L).
    """

    print("="*80)
    print("PARTIAL STRUCTURE ANALYSIS")
    print("="*80)
    print()

    # Smaller test: K=1, 4×4 matrix
    print("Testing K=1 (4×4 matrix) for tractability...")
    print()

    z, L, N = symbols('z L N', complex=True)
    z_conj = conjugate(z)

    def v_plus(n_val, z, L, N):
        freq_term = 2*pi*(n_val - N) / L
        numer_val = sqrt(2) * sin(L*z / 2)
        denom_val = L * (z + freq_term)
        return numer_val / denom_val

    def v_minus(n_val, z, L, N):
        freq_term = 2*pi*(n_val - N) / L
        numer_val = sqrt(2) * sin(L*z / 2)
        denom_val = L * (z - freq_term)
        return numer_val / denom_val

    matrix_rows = []
    for n_val in range(4):
        row = [v_plus(n_val, z, L, N), v_minus(n_val, z, L, N),
               v_plus(n_val, z_conj, L, N), v_minus(n_val, z_conj, L, N)]
        matrix_rows.append(row)

    M = Matrix(matrix_rows)

    print("K=1 (4×4) matrix structure:")
    print(f"  Shape: {M.shape}")
    print()

    print("Computing K=1 determinant...")
    sys.stdout.flush()

    try:
        d = det(M)
        print(f"K=1 determinant computed. Size: {len(str(d))} chars")
        print()

        # Try trigsimp
        print("Applying trigsimp (trig simplification)...")
        d_trig = trigsimp(d)
        print(f"After trigsimp: {len(str(d_trig))} chars")
        print()

        # Try to factor
        print("Attempting factorization...")
        d_factored = factor(d_trig)
        print(f"Factored: {len(str(d_factored))} chars")
        print()
        print("K=1 determinant (factored):")
        print(d_factored)
        print()

        return d_factored

    except Exception as e:
        print(f"ERROR: {e}")
        return None


if __name__ == '__main__':
    # Try full K=2 first
    result, syms = attempt_symbolic_factorization()

    if result is None:
        # Fall back to K=1
        print()
        print("Falling back to K=1 (smaller) case for structure analysis...")
        print()
        result = analyze_partial_structure()

    print()
    print("="*80)
    print("END SYMBOLIC ANALYSIS")
    print("="*80)

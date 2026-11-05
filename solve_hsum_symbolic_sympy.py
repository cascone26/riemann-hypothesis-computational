"""
EXACT SYMBOLIC APPROACH (sympy):

Build hsum constraint matrix with EXACT entries (Rational, not float).
Use sympy's exact nullspace() to find a basis vector.
Evaluate cubic/quartic polynomial exactly with sympy.simplify().

This avoids all floating-point precision issues completely.
"""

from sympy import (
    symbols, sqrt, sin, cos, I, Rational, Matrix, legendre,
    simplify, nsimplify, expand, summation, latex
)
from sympy import pi as sym_pi
import sympy as sp

def vPlusLC_symbolic(n, L, z):
    """vPlus_{LC}(n) = sqrt(2)*sin(L*z/2) * P_n(z) symbolically"""
    sqrt2 = sqrt(2)
    sin_term = sin(L * z / 2)
    leg_n = legendre(n, z)  # Legendre polynomial
    return sqrt2 * sin_term * leg_n

def vMinusLC_symbolic(n, L, z):
    """vMinus_{LC}(n) = sqrt(2)*sin(L*z/2) * U_n(z) symbolically"""
    sqrt2 = sqrt(2)
    sin_term = sin(L * z / 2)
    # U_n(z) is Chebyshev polynomial of the second kind
    from sympy import chebyshevu
    cheb_n = chebyshevu(n, z)
    return sqrt2 * sin_term * cheb_n

def build_hsum_matrix_symbolic(K, z_values, L, num_basis=10):
    """Build hsum matrix with exact symbolic entries."""
    num_rows = num_basis + 1
    num_cols = 4 * K

    M_list = []

    for n in range(num_basis + 1):
        row = []
        for i in range(K):
            z_i = z_values[i]
            z_i_conj = sp.conjugate(z_i)

            v_plus_i = vPlusLC_symbolic(n, L, z_i)
            v_minus_i = vMinusLC_symbolic(n, L, z_i)
            v_plus_conj_i = vPlusLC_symbolic(n, L, z_i_conj)
            v_minus_conj_i = vMinusLC_symbolic(n, L, z_i_conj)

            row.extend([v_plus_i, v_minus_i, v_plus_conj_i, v_minus_conj_i])

        M_list.append(row)

    M = Matrix(M_list)
    return M

def find_nullspace_symbolic(M):
    """Find exact nullspace basis using sympy.nullspace()."""
    print("  Computing nullspace (exact symbolic)...")
    NS = M.nullspace()
    print(f"  Nullspace dimension: {len(NS)}")

    if len(NS) == 0:
        return None

    # Return first basis vector
    return NS[0]

def verify_solution_symbolic(M, x):
    """Verify: M @ x = 0 symbolically."""
    residual = M * x
    print(f"  Residual vector (M @ x):")

    simplified_residuals = []
    for i in range(residual.rows):
        res = residual[i, 0]
        print(f"    [{i}] = {res}")

        # Try to simplify
        res_simplified = simplify(res)
        simplified_residuals.append(res_simplified)

        if res_simplified != res:
            print(f"           simplified = {res_simplified}")

    # Check if all simplified entries are zero
    all_zero = all(r == 0 for r in simplified_residuals)
    return all_zero

if __name__ == "__main__":
    print("=" * 120)
    print("K=5 CONTROL: Exact symbolic computation with sympy")
    print("=" * 120)

    # Simplified family: [0, 1, 2, 3, 4] + (5/2)j (all exact Rational)
    # This is cleaner than Powers-of-2 for symbolic computation
    z_test_symbolic = [
        0 + Rational(5, 2) * I,
        1 + Rational(5, 2) * I,
        2 + Rational(5, 2) * I,
        3 + Rational(5, 2) * I,
        4 + Rational(5, 2) * I
    ]

    K_control = 5
    L = 20  # Exact integer
    num_basis = 2  # Very small for manageable symbolic computation

    print(f"\nFamily: [0, 1, 2, 3, 4] + (5/2)j (exact Rational)")
    print(f"K = {K_control}, L = {L}, num_basis = {num_basis}")

    # Build hsum matrix
    print(f"\nBuilding hsum matrix with exact symbolic entries...")
    M_k5 = build_hsum_matrix_symbolic(K_control, z_test_symbolic, L, num_basis=num_basis)
    print(f"Matrix shape: {M_k5.shape}")
    print(f"Matrix entries (sample): {M_k5[0, 0]}")

    # Find nullspace
    print(f"\nFinding exact nullspace...")
    try:
        x_k5 = find_nullspace_symbolic(M_k5)

        if x_k5 is None:
            print("ERROR: No nullspace found")
        else:
            print(f"Nullspace solution found:")
            print(f"  Vector dimension: {len(x_k5)}")
            print(f"  First 3 entries: {x_k5[:3]}")

            # Verify
            print(f"\nVerifying solution...")
            is_valid = verify_solution_symbolic(M_k5, x_k5)

            if is_valid:
                print(f"✓ Solution is EXACT and CORRECT (M @ x = 0)")
            else:
                print(f"✗ Solution is WRONG (residuals not zero)")

            # Extract coefficients
            print(f"\nCoefficients from solution:")
            print(f"  p_0 (x[0]) = {x_k5[0]}")
            print(f"  q_0 (x[1]) = {x_k5[1]}")
            print(f"  r_0 (x[2]) = {x_k5[2]}")
            print(f"  s_0 (x[3]) = {x_k5[3]}")

    except Exception as e:
        print(f"ERROR during nullspace computation: {e}")
        import traceback
        traceback.print_exc()

    print(f"\n{'='*120}")

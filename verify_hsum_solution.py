"""
Debug: Verify that the extracted nullspace solution actually satisfies hsum=0.

If a solution extracted from nullspace doesn't satisfy the constraint,
there's a bug in how we're building or extracting the nullspace.
"""

import mpmath as mp
import numpy as np
from scipy.linalg import svd

def vPlusLC(n, L, z):
    if isinstance(n, int):
        n = mp.mpf(n)
    if isinstance(L, float):
        L = mp.mpf(L)
    if isinstance(z, complex):
        z = mp.mpc(z)
    sqrt2 = mp.sqrt(mp.mpf(2))
    arg = L * z / mp.mpf(2)
    sin_term = mp.sin(arg)
    leg_n = mp.legendre(n, z)
    return sqrt2 * sin_term * leg_n

def vMinusLC(n, L, z):
    if isinstance(n, int):
        n = mp.mpf(n)
    if isinstance(L, float):
        L = mp.mpf(L)
    if isinstance(z, complex):
        z = mp.mpc(z)
    sqrt2 = mp.sqrt(mp.mpf(2))
    arg = L * z / mp.mpf(2)
    sin_term = mp.sin(arg)
    cheb_n = mp.chebyu(n, z)
    return sqrt2 * sin_term * cheb_n

def build_hsum_matrix(K, z_values, L, num_basis=10):
    num_cols = 4 * K
    num_rows = num_basis + 1

    matrix = []
    for n in range(num_basis + 1):
        row = []
        for i in range(K):
            z_i = z_values[i] if i < len(z_values) else mp.mpc(0)
            z_i_conj = mp.conj(z_i)

            v_plus_i = vPlusLC(n, L, z_i)
            v_minus_i = vMinusLC(n, L, z_i)
            v_plus_conj_i = vPlusLC(n, L, z_i_conj)
            v_minus_conj_i = vMinusLC(n, L, z_i_conj)

            row.extend([v_plus_i, v_minus_i, v_plus_conj_i, v_minus_conj_i])

        matrix.append(row)

    # Convert to float64
    matrix_float = np.array([[complex(entry) for entry in row]
                             for row in matrix], dtype=np.complex128)
    return matrix_float, matrix  # Return both float and mpmath versions

def verify_solution(matrix_mpmath, solution_mpmath):
    """Check: matrix @ solution = 0?"""
    K = len(solution_mpmath) // 4
    num_basis = len(matrix_mpmath) - 1

    residuals = []
    for n in range(num_basis + 1):
        row = matrix_mpmath[n]
        dot_product = mp.mpc(0)
        for j, coeff in enumerate(solution_mpmath):
            dot_product += row[j] * coeff
        residuals.append(dot_product)

    return residuals  # Should all be ~0 if solution is correct

if __name__ == "__main__":
    mp.dps = 50
    L = mp.mpf(20)

    # K=5 control test
    z_powers_of_2 = [mp.mpc(float(2**k), 2.3) for k in range(6)]
    z_control = z_powers_of_2[:5]

    K_control = 5

    print("=" * 120)
    print("VERIFICATION: Does extracted solution satisfy hsum=0?")
    print("=" * 120)
    print(f"\nK=5 control test on Powers-of-2 family")

    # Build matrix
    matrix_float, matrix_mpmath = build_hsum_matrix(K_control, z_control, L, num_basis=10)

    # SVD and nullspace
    U, s, Vh = svd(matrix_float, full_matrices=True)
    rank = np.sum(s > 1e-10)
    nullspace_dim = matrix_float.shape[1] - rank

    print(f"\nhsum matrix shape: {matrix_float.shape}")
    print(f"Rank: {rank}, Nullspace dimension: {nullspace_dim}")

    # Extract first nullspace vector
    nullspace_basis = Vh[rank:, :].conj().T  # (4K, nullspace_dim)
    solution_float = nullspace_basis[:, 0]

    # Convert to mpmath
    solution_mpmath = np.array([mp.mpc(c) for c in solution_float])

    print(f"\nExtracted solution vector (first 5 entries):")
    print(f"  {[mp.nstr(s, 10) for s in solution_mpmath[:5]]}")

    # Verify: does matrix @ solution = 0?
    residuals = verify_solution(matrix_mpmath, solution_mpmath)

    print(f"\nResiduals (matrix @ solution for each constraint n=0..10):")
    for i, res in enumerate(residuals):
        print(f"  n={i}: {mp.nstr(mp.fabs(res), 10)}")

    max_residual = max(mp.fabs(r) for r in residuals)
    print(f"\nMax residual: {mp.nstr(max_residual, 15)}")

    if max_residual < mp.mpf(1e-30):
        print("\n✓ SOLUTION IS CORRECT: Satisfies hsum=0 at high precision")
    elif max_residual < mp.mpf(1e-10):
        print("\n⚠ SOLUTION IS MARGINAL: Satisfies hsum≈0 at float precision only")
    else:
        print("\n✗ SOLUTION IS WRONG: Does NOT satisfy hsum (debug needed)")

    print(f"\n{'='*120}")

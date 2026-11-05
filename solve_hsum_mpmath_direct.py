"""
FIX: Use mpmath's own linear algebra (not scipy) to solve hsum constraint.

scipy's SVD fails on large-scale ill-conditioned matrices (float64 inadequate).
mpmath can handle arbitrary precision and will give us reliable nullspace.
"""

import mpmath as mp
import numpy as np

def vPlusLC(n, L, z):
    if isinstance(n, int):
        n = mp.mpf(n)
    if isinstance(L, (float, int)):
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
    if isinstance(L, (float, int)):
        L = mp.mpf(L)
    if isinstance(z, complex):
        z = mp.mpc(z)
    sqrt2 = mp.sqrt(mp.mpf(2))
    arg = L * z / mp.mpf(2)
    sin_term = mp.sin(arg)
    cheb_n = mp.chebyu(n, z)
    return sqrt2 * sin_term * cheb_n

def build_hsum_matrix_mpmath(K, z_values, L, num_basis=10):
    """Build hsum matrix using mpmath (arbitrary precision)."""
    num_rows = num_basis + 1
    num_cols = 4 * K

    # Use mpmath matrix
    M = mp.matrix(num_rows, num_cols)

    for n in range(num_basis + 1):
        col_idx = 0
        for i in range(K):
            z_i = z_values[i] if i < len(z_values) else mp.mpc(0)
            z_i_conj = mp.conj(z_i)

            v_plus_i = vPlusLC(n, L, z_i)
            v_minus_i = vMinusLC(n, L, z_i)
            v_plus_conj_i = vPlusLC(n, L, z_i_conj)
            v_minus_conj_i = vMinusLC(n, L, z_i_conj)

            M[n, col_idx] = v_plus_i
            M[n, col_idx + 1] = v_minus_i
            M[n, col_idx + 2] = v_plus_conj_i
            M[n, col_idx + 3] = v_minus_conj_i

            col_idx += 4

    return M

def find_nullspace_mpmath(M):
    """Find one vector in the nullspace of M using mpmath LU decomposition."""
    # Strategy: Set one variable to 1, solve for others
    # This is simpler than full nullspace computation but gives a solution

    # Try setting the last variable to 1 and solving for the rest
    num_rows, num_cols = M.rows, M.cols

    # Create augmented system: M_reduced @ x_reduced = -M[:, -1] * 1
    M_reduced = mp.matrix(num_rows, num_cols - 1)
    b = mp.matrix(num_rows, 1)

    for i in range(num_rows):
        for j in range(num_cols - 1):
            M_reduced[i, j] = M[i, j]
        b[i, 0] = -M[i, num_cols - 1]  # Last column, coefficient = 1

    try:
        # Solve M_reduced @ x_reduced = b
        x_reduced = mp.lu_solve(M_reduced, b)

        # Construct full solution: [x_reduced; 1]
        solution = mp.matrix(num_cols, 1)
        for i in range(num_cols - 1):
            solution[i, 0] = x_reduced[i, 0]
        solution[num_cols - 1, 0] = mp.mpc(1)

        return solution
    except:
        # If that fails, try a different approach: set first variable to 1
        M_reduced = mp.matrix(num_rows, num_cols - 1)
        b = mp.matrix(num_rows, 1)

        for i in range(num_rows):
            for j in range(1, num_cols):
                M_reduced[i, j - 1] = M[i, j]
            b[i, 0] = -M[i, 0]  # First column, coefficient = 1

        try:
            x_reduced = mp.lu_solve(M_reduced, b)
            solution = mp.matrix(num_cols, 1)
            solution[0, 0] = mp.mpc(1)
            for i in range(num_cols - 1):
                solution[i + 1, 0] = x_reduced[i, 0]
            return solution
        except:
            return None

def verify_solution(M, x):
    """Check: M @ x = 0?"""
    residual = M * x
    max_res = max(mp.fabs(residual[i, 0]) for i in range(residual.rows))
    return residual, max_res

if __name__ == "__main__":
    mp.dps = 50
    L = mp.mpf(20)

    # K=5 control test
    z_powers_of_2 = [mp.mpc(float(2**k), 2.3) for k in range(6)]
    z_control = z_powers_of_2[:5]

    K_control = 5

    print("=" * 120)
    print("K=5 CONTROL TEST: mpmath nullspace solver")
    print("=" * 120)

    # Build matrix with mpmath
    print(f"\nBuilding hsum matrix for K=5...")
    M_k5 = build_hsum_matrix_mpmath(K_control, z_control, L, num_basis=10)
    print(f"Matrix shape: {M_k5.rows} × {M_k5.cols}")

    # Find solution
    print(f"Finding nullspace solution...")
    x_k5 = find_nullspace_mpmath(M_k5)

    if x_k5 is None:
        print("ERROR: Could not find solution")
    else:
        print(f"Solution found (first 5 entries):")
        for i in range(min(5, x_k5.rows)):
            print(f"  x[{i}] = {mp.nstr(x_k5[i, 0], 15)}")

        # Verify
        print(f"\nVerifying: M @ x = 0?")
        residual, max_res = verify_solution(M_k5, x_k5)

        print(f"Residuals for each constraint n=0..10:")
        for i in range(residual.rows):
            print(f"  n={i}: {mp.nstr(mp.fabs(residual[i, 0]), 12)}")

        print(f"\nMax residual: {mp.nstr(max_res, 15)}")

        if max_res < mp.mpf(1e-40):
            print("✓ Solution is CORRECT (satisfies hsum=0)")
        else:
            print("✗ Solution is WRONG (does not satisfy hsum)")

    print(f"\n{'='*120}")

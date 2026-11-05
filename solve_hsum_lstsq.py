"""
ANOTHER APPROACH: Use scipy.linalg.lstsq (least-squares) which is more robust.

lstsq finds the minimum-norm solution to the underdetermined system.
It should handle ill-conditioning better than manual SVD.
"""

import mpmath as mp
import numpy as np
from scipy.linalg import lstsq

def vPlusLC(n, L, z, dps):
    """vPlus_{LC}(n) with explicit dps setting"""
    mp.dps = dps
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

def vMinusLC(n, L, z, dps):
    """vMinus_{LC}(n) with explicit dps setting"""
    mp.dps = dps
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

def build_hsum_matrix_float(K, z_values, L, num_basis=10):
    """Build hsum matrix with high precision mpmath, then convert to float64."""
    dps = 100
    mp.dps = dps

    num_rows = num_basis + 1
    num_cols = 4 * K

    matrix_list = []

    for n in range(num_basis + 1):
        row = []
        for i in range(K):
            z_i = z_values[i] if i < len(z_values) else mp.mpc(0)
            z_i_conj = mp.conj(z_i)

            v_plus_i = vPlusLC(n, L, z_i, dps)
            v_minus_i = vMinusLC(n, L, z_i, dps)
            v_plus_conj_i = vPlusLC(n, L, z_i_conj, dps)
            v_minus_conj_i = vMinusLC(n, L, z_i_conj, dps)

            # Convert to complex float (accepting precision loss)
            row.extend([
                complex(v_plus_i),
                complex(v_minus_i),
                complex(v_plus_conj_i),
                complex(v_minus_conj_i)
            ])

        matrix_list.append(row)

    matrix_float = np.array(matrix_list, dtype=np.complex128)
    return matrix_float

def find_solution_lstsq(M):
    """Use lstsq to find minimum-norm solution to M @ x = 0 (underdetermined)."""
    # For homogeneous system M @ x = 0, lstsq will find the minimum-norm solution
    # This is equivalent to finding the last singular vector (corresponding to zero singular value)

    b = np.zeros(M.shape[0], dtype=np.complex128)  # Right-hand side = 0
    x, residuals, rank, s = lstsq(M, b, lapack_driver='gelsy')

    print(f"  lstsq results:")
    print(f"    Rank of M: {rank}")
    print(f"    Residuals (if any): {residuals}")
    if s is not None:
        print(f"    Singular values (first 5): {s[:5] if len(s) >= 5 else s}")
    else:
        print(f"    Singular values: (not returned)")

    return x

def verify_solution_mpmath(M_mpmath, x_float):
    """Verify solution using mpmath (high precision)."""
    # Convert x to mpmath
    x_mp = np.array([mp.mpc(c) for c in x_float])

    # Compute M @ x
    residual = M_mpmath * mp.matrix(x_mp.tolist())

    residuals_list = [residual[i, 0] for i in range(residual.rows)]
    max_res = max(mp.fabs(r) for r in residuals_list)

    return residuals_list, max_res

def build_hsum_matrix_mpmath(K, z_values, L, num_basis=10):
    """Build hsum matrix using mpmath for verification."""
    dps = 100
    mp.dps = dps

    num_rows = num_basis + 1
    num_cols = 4 * K

    M = mp.matrix(num_rows, num_cols)

    for n in range(num_basis + 1):
        col_idx = 0
        for i in range(K):
            z_i = z_values[i] if i < len(z_values) else mp.mpc(0)
            z_i_conj = mp.conj(z_i)

            v_plus_i = vPlusLC(n, L, z_i, dps)
            v_minus_i = vMinusLC(n, L, z_i, dps)
            v_plus_conj_i = vPlusLC(n, L, z_i_conj, dps)
            v_minus_conj_i = vMinusLC(n, L, z_i_conj, dps)

            M[n, col_idx] = v_plus_i
            M[n, col_idx + 1] = v_minus_i
            M[n, col_idx + 2] = v_plus_conj_i
            M[n, col_idx + 3] = v_minus_conj_i

            col_idx += 4

    return M

if __name__ == "__main__":
    # K=5 control test
    z_powers_of_2 = [mp.mpc(float(2**k), 2.3) for k in range(6)]
    z_control = z_powers_of_2[:5]

    K_control = 5
    L = mp.mpf(20)

    print("=" * 120)
    print("K=5 CONTROL: Using lstsq to find hsum=0 solution")
    print("=" * 120)

    # Build matrix with high precision mpmath, convert to float64
    print(f"\nBuilding hsum matrix (dps=100)...")
    M_float = build_hsum_matrix_float(K_control, z_control, L, num_basis=10)
    print(f"Matrix shape: {M_float.shape}")
    print(f"Condition number: {np.linalg.cond(M_float):.2e}")

    # Find solution using lstsq
    print(f"\nSolving M @ x = 0 using lstsq...")
    x_float = find_solution_lstsq(M_float)

    print(f"\nSolution vector (first 5 entries):")
    for i in range(min(5, len(x_float))):
        print(f"  x[{i}] = {x_float[i]:.10e}")

    # Verify using mpmath
    print(f"\nVerifying solution using mpmath (dps=100)...")
    M_mpmath = build_hsum_matrix_mpmath(K_control, z_control, L, num_basis=10)
    residuals, max_res = verify_solution_mpmath(M_mpmath, x_float)

    print(f"Residuals for each constraint n=0..10:")
    for i, res in enumerate(residuals):
        print(f"  n={i}: {mp.nstr(mp.fabs(res), 12)}")

    print(f"\nMax residual: {mp.nstr(max_res, 15)}")

    if max_res < mp.mpf(1e-30):
        print("✓ Solution is CORRECT")
    else:
        print("✗ Solution still doesn't satisfy hsum (precision/method issue)")

    print(f"\n{'='*120}")

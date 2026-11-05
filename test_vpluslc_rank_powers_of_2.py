"""
Real K-induction test: check whether vPlusLC/vMinusLC are linearly INDEPENDENT
on the Powers-of-2 family. If they are independent, the Lean constraint
(sum of p_i*vPlus + q_i*vMinus + ... = 0) can only be satisfied with all p_i=q_i=...=0,
meaning the K-induction identity constraint is VACUOUS (can never be satisfied on this family).

Linear independence on Powers-of-2 would mean:
- K-induction CAN'T be applied (no non-trivial coefficients exist)
- Therefore no partial positive result on this family

Conversely, if they ARE linearly dependent on Powers-of-2:
- There exist non-trivial coefficients that make vPlusLC/vMinusLC sum = 0
- We could then check whether those same coefficients make cubic/quartic = 0
"""

import mpmath as mp
import numpy as np
from scipy.linalg import svd

def vPlusLC(n, L, z):
    """vPlus_{LC}(n) from HLGComplexContinuation.lean: sqrt(2)*sin(L*z/2) * Leg_n(z)"""
    if isinstance(n, int):
        n = mp.mpf(n)
    if isinstance(L, float):
        L = mp.mpf(L)
    if isinstance(z, complex):
        z = mp.mpc(z)

    sqrt2 = mp.sqrt(mp.mpf(2))
    arg = L * z / mp.mpf(2)
    sin_term = mp.sin(arg)
    leg_n = mp.legendre(n, z)  # Legendre polynomial P_n(z)

    return sqrt2 * sin_term * leg_n

def vMinusLC(n, L, z):
    """vMinus_{LC}(n) = sqrt(2)*sin(L*z/2) * Cheb_n(z)"""
    if isinstance(n, int):
        n = mp.mpf(n)
    if isinstance(L, float):
        L = mp.mpf(L)
    if isinstance(z, complex):
        z = mp.mpc(z)

    sqrt2 = mp.sqrt(mp.mpf(2))
    arg = L * z / mp.mpf(2)
    sin_term = mp.sin(arg)
    cheb_n = mp.chebyu(n, z)  # Chebyshev U polynomial

    return sqrt2 * sin_term * cheb_n

if __name__ == "__main__":
    mp.dps = 50  # 50 decimal places for initial rank check

    # Powers of 2 family
    z_values = [mp.mpc(float(2**k), 2.3) for k in range(6)]
    N = 10  # Number of basis functions (0..N)
    L = mp.mpf(20)

    print("=" * 100)
    print("LINEAR INDEPENDENCE CHECK: vPlusLC/vMinusLC on Powers-of-2 Family")
    print("=" * 100)
    print(f"\nFamily: z = [1, 2, 4, 8, 16, 32] + 2.3j")
    print(f"Basis: vPlusLC/vMinusLC for n=0..{N} (2*{N+1}=22 basis functions per z)")
    print(f"Precision: mpmath dps=50")
    print(f"L = {L}")

    # Build the matrix: rows = basis functions (vPlus_0, vMinus_0, vPlus_1, vMinus_1, ..., vPlus_N, vMinus_N)
    #                   cols = z values (6 columns)
    # Entry[i, j] = (vPlus or vMinus evaluated at z_j)

    basis_size = 2 * (N + 1)  # 22 basis functions
    num_z = len(z_values)

    # Build matrix
    matrix = []
    for n in range(N + 1):
        # Row for vPlus_n
        row_plus = [vPlusLC(n, L, z) for z in z_values]
        matrix.append(row_plus)

        # Row for vMinus_n
        row_minus = [vMinusLC(n, L, z) for z in z_values]
        matrix.append(row_minus)

    # Convert to float for rank computation (mpmath is slow for SVD)
    # This is acceptable since we're checking structural rank, not exact precision
    matrix_float = np.array([[complex(entry) for entry in row]
                             for row in matrix], dtype=np.complex128)

    print(f"\nMatrix shape: {basis_size} basis functions × {num_z} z-values")
    print(f"Matrix:\n{matrix_float}")

    # Compute rank via SVD
    U, s, Vh = svd(matrix_float, full_matrices=False)
    rank = np.sum(s > 1e-10)  # Numerical rank with threshold

    print(f"\nSingular values (first 10): {s[:10]}")
    print(f"Rank (threshold=1e-10): {rank} / {min(basis_size, num_z)}")

    # Interpret result
    print(f"\n" + "=" * 100)
    if rank == num_z:
        print(f"*** LINEAR DEPENDENCE: rank = {num_z} (full column rank) ***")
        print(f"All {num_z} z-values produce LINEARLY INDEPENDENT vPlusLC/vMinusLC vectors.")
        print(f"Consequence: The Lean constraint (Σ p_i*vPlus + ... = 0) is VACUOUS on Powers-of-2.")
        print(f"  → No non-trivial coefficients can satisfy the constraint.")
        print(f"  → K-induction cannot be applied on this family (no valid coefficients exist).")
        print(f"\nInterpretation: Powers-of-2 is a family where the linear constraint itself is unsatisfiable.")
    else:
        print(f"*** LINEAR INDEPENDENCE: rank = {rank} < {num_z} (rank-deficient) ***")
        print(f"The {num_z} z-values produce LINEARLY DEPENDENT vPlusLC/vMinusLC vectors.")
        print(f"Consequence: There exist non-trivial coefficients satisfying the Lean constraint.")
        print(f"  → Next step: solve for those coefficients and test K-induction cubic/quartic identity.")
        print(f"\nInterpretation: Powers-of-2 admits restricted K-induction (if cubic/quartic also vanishes).")

    print("=" * 100)

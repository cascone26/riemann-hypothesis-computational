"""
Rank check for π-multiples family: z = [π, 2π, 3π, 4π, 5π, 6π] + 3j
Same test as Powers-of-2 but on a different survivor family.
"""

import mpmath as mp
import numpy as np
from scipy.linalg import svd

def vPlusLC(n, L, z):
    """vPlus_{LC}(n) = sqrt(2)*sin(L*z/2) * Leg_n(z)"""
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
    cheb_n = mp.chebyu(n, z)

    return sqrt2 * sin_term * cheb_n

if __name__ == "__main__":
    mp.dps = 50

    # π-multiples family
    z_values = [mp.mpc(float(k * np.pi), 3.0) for k in range(1, 7)]
    N = 10
    L = mp.mpf(20)

    print("=" * 100)
    print("LINEAR INDEPENDENCE CHECK: vPlusLC/vMinusLC on π-multiples Family")
    print("=" * 100)
    print(f"\nFamily: z = [π, 2π, 3π, 4π, 5π, 6π] + 3j")
    print(f"Precision: mpmath dps=50")
    print(f"L = {L}")

    # Build matrix
    basis_size = 2 * (N + 1)
    num_z = len(z_values)

    matrix = []
    for n in range(N + 1):
        row_plus = [vPlusLC(n, L, z) for z in z_values]
        matrix.append(row_plus)

        row_minus = [vMinusLC(n, L, z) for z in z_values]
        matrix.append(row_minus)

    # Convert to float
    matrix_float = np.array([[complex(entry) for entry in row]
                             for row in matrix], dtype=np.complex128)

    print(f"\nMatrix shape: {basis_size} basis functions × {num_z} z-values")

    # Compute rank
    U, s, Vh = svd(matrix_float, full_matrices=False)
    rank = np.sum(s > 1e-10)

    print(f"Singular values (first 10): {s[:10]}")
    print(f"Rank (threshold=1e-10): {rank} / {min(basis_size, num_z)}")

    # Interpret
    print(f"\n" + "=" * 100)
    if rank == num_z:
        print(f"*** LINEAR INDEPENDENCE: rank = {num_z} (full column rank) ***")
        print(f"Result: π-multiples ALSO has linearly independent vPlusLC/vMinusLC vectors.")
        print(f"Conclusion: The constraint obstruction is UNIVERSAL across survivors so far.")
    elif rank < num_z:
        print(f"*** LINEAR DEPENDENCE: rank = {rank} < {num_z} (rank-deficient) ***")
        print(f"Result: π-multiples ADMITS non-trivial constraint coefficients.")
        print(f"Next: Test whether cubic/quartic sum vanishes with those coefficients.")

    print("=" * 100)

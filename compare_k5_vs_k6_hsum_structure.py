"""
STRUCTURAL COMPARISON: hsum constraint matrix for K=5 (known-working) vs K=6 (survivors).

K=5's working constructions have independent vPlusLC/vMinusLC vectors (S2693).
Yet K=5 doesn't show vacuous-only constraints. Why would K=6 be structurally different?

Hypothesis: The hsum constraint's STRUCTURE (not vPlusLC/vMinusLC independence) differs between K=5 and K=6,
making the solution space vanish at K=6 even when vectors are independent.

Test: For both K=5 and K=6, compute the hsum constraint matrix and check:
- Rank (solution-space dimension)
- Nullspace structure (what solutions exist?)

Using the SAME family (z-values) for both K to isolate the K-dependent structural difference.
"""

import mpmath as mp
import numpy as np
from scipy.linalg import svd, null_space

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

def build_hsum_matrix(K, z_values, L, num_basis=None):
    """
    Build the hsum constraint matrix.

    For each z_i in the family, we have 4 unknowns: p_i, q_i, r_i, s_i
    Total unknowns: 4K

    Constraints (rows): For each basis function n (from 0 to num_basis),
      Σ_i (p_i*vPlus(z_i) + q_i*vMinus(z_i) + r_i*vPlus(conj z_i) + s_i*vMinus(conj z_i)) = 0
    Total constraints: (num_basis + 1)

    Matrix shape: (num_basis + 1) rows × 4K columns
    Each row: [p_0, q_0, r_0, s_0, p_1, q_1, r_1, s_1, ...]
    Entry[n, 4i+0] = vPlus(n, z_i)
    Entry[n, 4i+1] = vMinus(n, z_i)
    Entry[n, 4i+2] = vPlus(n, conj z_i)
    Entry[n, 4i+3] = vMinus(n, conj z_i)
    """
    if num_basis is None:
        num_basis = 10

    num_z = len(z_values)
    num_cols = 4 * K  # p, q, r, s for each z_value
    num_rows = num_basis + 1  # n = 0 to num_basis

    matrix = []
    for n in range(num_basis + 1):
        row = []
        for i in range(K):
            z_i = z_values[i] if i < num_z else mp.mpc(0)  # Pad if needed
            z_i_conj = mp.conj(z_i)

            v_plus_i = vPlusLC(n, L, z_i)
            v_minus_i = vMinusLC(n, L, z_i)
            v_plus_conj_i = vPlusLC(n, L, z_i_conj)
            v_minus_conj_i = vMinusLC(n, L, z_i_conj)

            row.extend([v_plus_i, v_minus_i, v_plus_conj_i, v_minus_conj_i])

        matrix.append(row)

    # Convert to float
    matrix_float = np.array([[complex(entry) for entry in row]
                             for row in matrix], dtype=np.complex128)

    return matrix_float

if __name__ == "__main__":
    mp.dps = 50
    L = mp.mpf(20)

    # Test family: Powers-of-2 (used in fork's tests)
    z_test = [mp.mpc(float(2**k), 2.3) for k in range(6)]

    print("=" * 120)
    print("STRUCTURAL COMPARISON: K=5 (KNOWN-WORKING) vs K=6 (SURVIVOR FAMILY)")
    print("=" * 120)
    print(f"\nTest family: Powers-of-2, z = [1, 2, 4, 8, 16, 32] + 2.3j")
    print(f"L = {L}, num_basis = 10")

    # K=5 case (known working)
    print(f"\n{'='*120}")
    print("K=5 (known-working construction)")
    print(f"{'='*120}")
    matrix_k5 = build_hsum_matrix(5, z_test, L, num_basis=10)
    print(f"hsum matrix shape: {matrix_k5.shape} (rows × cols = constraints × unknowns)")
    print(f"  Constraints: 11 (basis n=0..10)")
    print(f"  Unknowns: 4*5 = 20 (p, q, r, s for each of 5 z-values)")

    U5, s5, Vh5 = svd(matrix_k5, full_matrices=False)
    rank_k5 = np.sum(s5 > 1e-10)
    nullspace_dim_k5 = 20 - rank_k5  # number of free parameters

    print(f"\nSingular values (first 10): {s5[:10]}")
    print(f"Rank: {rank_k5} / {min(matrix_k5.shape)}")
    print(f"Nullspace dimension: {nullspace_dim_k5}")
    print(f"  → Solution space is {nullspace_dim_k5}-dimensional (non-trivial solutions exist)")

    # K=6 case (fork's survivor)
    print(f"\n{'='*120}")
    print("K=6 (survivor family: Powers-of-2)")
    print(f"{'='*120}")
    matrix_k6 = build_hsum_matrix(6, z_test, L, num_basis=10)
    print(f"hsum matrix shape: {matrix_k6.shape} (rows × cols = constraints × unknowns)")
    print(f"  Constraints: 11 (basis n=0..10)")
    print(f"  Unknowns: 4*6 = 24 (p, q, r, s for each of 6 z-values)")

    U6, s6, Vh6 = svd(matrix_k6, full_matrices=False)
    rank_k6 = np.sum(s6 > 1e-10)
    nullspace_dim_k6 = 24 - rank_k6

    print(f"\nSingular values (first 10): {s6[:10]}")
    print(f"Rank: {rank_k6} / {min(matrix_k6.shape)}")
    print(f"Nullspace dimension: {nullspace_dim_k6}")
    if nullspace_dim_k6 == 0:
        print(f"  → ONLY TRIVIAL SOLUTION (constraint system is FULL RANK)")
    else:
        print(f"  → Solution space is {nullspace_dim_k6}-dimensional (non-trivial solutions exist)")

    # Structural comparison
    print(f"\n{'='*120}")
    print("STRUCTURAL INTERPRETATION")
    print(f"{'='*120}")
    print(f"\nK=5: {rank_k5} constraints rank / {20} unknowns → nullspace dim = {nullspace_dim_k5}")
    print(f"K=6: {rank_k6} constraints rank / {24} unknowns → nullspace dim = {nullspace_dim_k6}")

    if nullspace_dim_k5 > 0 and nullspace_dim_k6 == 0:
        print(f"\n*** STRUCTURAL DIFFERENCE FOUND ***")
        print(f"K=5: Non-trivial solutions exist ({nullspace_dim_k5}-dimensional solution space)")
        print(f"K=6: Only trivial solution (full-rank constraint system)")
        print(f"\nThis is NOT due to vPlusLC/vMinusLC independence (both K=5/K=6 have it),")
        print(f"but due to the hsum constraint's STRUCTURE CHANGING at K=6.")
        print(f"Hypothesis: The constraint matrix has different rank properties by K-dependent design.")
    elif nullspace_dim_k5 == 0 and nullspace_dim_k6 == 0:
        print(f"\n*** BOTH K=5 and K=6 have ONLY TRIVIAL SOLUTIONS ***")
        print(f"If K=5 works despite this, it means the K=5 theorem doesn't actually NEED")
        print(f"non-trivial hsum solutions — the vacuous constraint is sufficient.")
    else:
        print(f"\n*** UNEXPECTED: Both systems have non-trivial solution spaces ***")
        print(f"If so, the fork's conclusion (K=6 vacuous, K=5 non-vacuous) is WRONG.")

    print(f"\n{'='*120}")

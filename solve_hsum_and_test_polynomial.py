"""
THE REAL K-INDUCTION TEST (Layer 4):

1. Build hsum constraint matrix for K (5 or 6)
2. Solve hsum = 0 explicitly → extract a non-trivial solution
3. Plug coefficients into the cubic/quartic polynomial sum
4. Evaluate at high precision (dps=50+)
5. Return: actual polynomial value

Control: K=5 should give ~0 (known working)
Test: K=6 on survivors should reveal if polynomial identity holds

Family: Powers-of-2 (cleanest survivor)
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

def build_hsum_matrix(K, z_values, L, num_basis=10):
    """Build hsum constraint matrix: (num_basis+1) × (4K) """
    num_z = len(z_values)
    num_cols = 4 * K
    num_rows = num_basis + 1

    matrix = []
    for n in range(num_basis + 1):
        row = []
        for i in range(K):
            z_i = z_values[i] if i < num_z else mp.mpc(0)
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
    return matrix_float

def extract_nullspace_solution(matrix_float):
    """Compute nullspace and extract one non-trivial solution vector."""
    # SVD to get nullspace
    U, s, Vh = svd(matrix_float, full_matrices=True)

    # Nullspace vectors are the right singular vectors corresponding to zero singular values
    rank = np.sum(s > 1e-10)
    nullspace_dim = matrix_float.shape[1] - rank

    if nullspace_dim == 0:
        return None, 0  # No non-trivial solution

    # Last nullspace_dim columns of Vh.H (or equivalently, last rows of Vh)
    nullspace_basis = Vh[rank:, :].conj().T  # Shape: (4K, nullspace_dim)

    # Extract first nullspace basis vector as the solution
    solution = nullspace_basis[:, 0]  # Shape: (4K,)

    return solution, nullspace_dim

def evaluate_cubic(K, n_val, N, L, z_values, p, q, r, s, A_vals, A_conj_vals):
    """Evaluate cubic sum: Σ_i cubic_i (with actual p, q, r, s coefficients)"""
    cubic_sum = mp.mpc(0)

    for i in range(K):
        z_i = z_values[i]
        z_i_conj = mp.conj(z_i)
        A_i = A_vals[i]
        A_i_conj = A_conj_vals[i]

        p_i = p[i]
        q_i = q[i]
        r_i = r[i]
        s_i = s[i]

        xi = mp.mpf(2) * mp.pi * (mp.mpf(n_val) - mp.mpf(N)) / mp.mpf(L)
        xi2 = xi ** 2
        xi3 = xi ** 3

        # cubic_i from CrossK6Key.lean
        term1 = (p_i + q_i) * A_i * z_i * (z_i_conj ** 2)
        term2 = (r_i + s_i) * A_i_conj * (z_i ** 2) * z_i_conj
        term3 = (q_i - p_i) * A_i * (z_i_conj ** 2) * xi
        term4 = (s_i - r_i) * A_i_conj * (z_i ** 2) * xi
        term5 = -(p_i + q_i) * A_i * z_i * xi2
        term6 = -(r_i + s_i) * A_i_conj * z_i_conj * xi2
        term7 = (p_i - q_i) * A_i * xi3
        term8 = (r_i - s_i) * A_i_conj * xi3

        cubic_i = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8
        cubic_sum += cubic_i

    return cubic_sum

def evaluate_quartic_product(K, n_val, N, L, z_values):
    """Evaluate quartic product: Π_i quartic_i"""
    quartic_product = mp.mpc(1)

    for i in range(K):
        z_i = z_values[i]
        z_i_conj = mp.conj(z_i)

        xi = mp.mpf(2) * mp.pi * (mp.mpf(n_val) - mp.mpf(N)) / mp.mpf(L)

        # quartic_i from CrossK6Key.lean
        quartic_i = (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)
        quartic_product *= quartic_i

    return quartic_product

def evaluate_cubic_quartic_sum(K, n_val, N, L, z_values, coeffs, A_vals, A_conj_vals):
    """
    Evaluate the full cubic/quartic sum from CrossK6Key.lean:
    cubic_1 * Q_2 * Q_3 * ... * Q_K + cubic_2 * Q_1 * Q_3 * ... * Q_K + ...

    where Q_i = quartic_i
    """
    # Extract coefficients
    p = [mp.mpc(coeffs[4*i+0]) for i in range(K)]
    q = [mp.mpc(coeffs[4*i+1]) for i in range(K)]
    r = [mp.mpc(coeffs[4*i+2]) for i in range(K)]
    s = [mp.mpc(coeffs[4*i+3]) for i in range(K)]

    # Compute all quartics
    quartics = []
    for i in range(K):
        z_i = z_values[i]
        z_i_conj = mp.conj(z_i)
        xi = mp.mpf(2) * mp.pi * (mp.mpf(n_val) - mp.mpf(N)) / mp.mpf(L)
        quartic_i = (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)
        quartics.append(quartic_i)

    # Compute all cubics
    cubics = []
    for i in range(K):
        z_i = z_values[i]
        z_i_conj = mp.conj(z_i)
        A_i = A_vals[i]
        A_i_conj = A_conj_vals[i]

        p_i = p[i]
        q_i = q[i]
        r_i = r[i]
        s_i = s[i]

        xi = mp.mpf(2) * mp.pi * (mp.mpf(n_val) - mp.mpf(N)) / mp.mpf(L)
        xi2 = xi ** 2
        xi3 = xi ** 3

        term1 = (p_i + q_i) * A_i * z_i * (z_i_conj ** 2)
        term2 = (r_i + s_i) * A_i_conj * (z_i ** 2) * z_i_conj
        term3 = (q_i - p_i) * A_i * (z_i_conj ** 2) * xi
        term4 = (s_i - r_i) * A_i_conj * (z_i ** 2) * xi
        term5 = -(p_i + q_i) * A_i * z_i * xi2
        term6 = -(r_i + s_i) * A_i_conj * z_i_conj * xi2
        term7 = (p_i - q_i) * A_i * xi3
        term8 = (r_i - s_i) * A_i_conj * xi3

        cubic_i = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8
        cubics.append(cubic_i)

    # Sum: cubic_i * (all quartics except quartic_i)
    total = mp.mpc(0)
    for i in range(K):
        product_others = mp.mpc(1)
        for j in range(K):
            if i != j:
                product_others *= quartics[j]
        total += cubics[i] * product_others

    return total

if __name__ == "__main__":
    mp.dps = 50
    L = mp.mpf(20)
    N = 10
    n_test = 15

    # Powers-of-2 family
    z_powers_of_2 = [mp.mpc(float(2**k), 2.3) for k in range(6)]

    print("=" * 120)
    print("LAYER 4: EXPLICIT hsum SOLUTION → CUBIC/QUARTIC EVALUATION")
    print("=" * 120)
    print(f"\nFamily: Powers-of-2, z = [1, 2, 4, 8, 16, 32] + 2.3j")
    print(f"L = {L}, N = {N}, test at n = {n_test}")

    # ========== CONTROL: K=5 ==========
    print(f"\n{'='*120}")
    print("CONTROL TEST: K=5 (known-working, should give ~0)")
    print(f"{'='*120}")

    K_control = 5
    z_control = z_powers_of_2[:K_control]  # Use first 5 z-values

    # Build hsum matrix for K=5
    matrix_k5 = build_hsum_matrix(K_control, z_control, L, num_basis=10)
    print(f"\nhsum matrix shape: {matrix_k5.shape}")
    print(f"Rank: {np.linalg.matrix_rank(matrix_k5)}")

    # Solve for non-trivial coefficients
    solution_k5_float, nullspace_dim_k5 = extract_nullspace_solution(matrix_k5)

    if solution_k5_float is None:
        print("ERROR: K=5 has no non-trivial solution (method bug)")
    else:
        print(f"Nullspace dimension: {nullspace_dim_k5}")
        print(f"Extracted solution vector (first 5 entries): {solution_k5_float[:5]}")

        # Convert to mpmath
        solution_k5 = np.array([mp.mpc(c) for c in solution_k5_float])

        # Compute A values for K=5
        sqrt2 = mp.sqrt(mp.mpf(2))
        A_k5 = [sqrt2 * mp.sin(L * z / mp.mpf(2)) for z in z_control]
        A_conj_k5 = [mp.conj(A) for A in A_k5]

        # Evaluate cubic/quartic sum for K=5
        result_k5 = evaluate_cubic_quartic_sum(K_control, n_test, N, L, z_control, solution_k5, A_k5, A_conj_k5)

        print(f"\nK=5 cubic/quartic sum with hsum solution:")
        print(f"  {mp.nstr(result_k5, 20)}")
        print(f"  |result| = {mp.nstr(mp.fabs(result_k5), 15)}")

        if mp.fabs(result_k5) < mp.mpf(1e-20):
            print(f"  ✓ K=5 WORKS: result ≈ 0 (method verified)")
            k5_works = True
        else:
            print(f"  ✗ K=5 FAILS: result ≠ 0 (method bug?)")
            k5_works = False

    # ========== TEST: K=6 ==========
    print(f"\n{'='*120}")
    print("TEST: K=6 (survivor family, what does the polynomial give?)")
    print(f"{'='*120}")

    K_test = 6

    # Build hsum matrix for K=6
    matrix_k6 = build_hsum_matrix(K_test, z_powers_of_2, L, num_basis=10)
    print(f"\nhsum matrix shape: {matrix_k6.shape}")
    print(f"Rank: {np.linalg.matrix_rank(matrix_k6)}")

    # Solve for non-trivial coefficients
    solution_k6_float, nullspace_dim_k6 = extract_nullspace_solution(matrix_k6)

    if solution_k6_float is None:
        print("ERROR: K=6 has no non-trivial solution")
    else:
        print(f"Nullspace dimension: {nullspace_dim_k6}")
        print(f"Extracted solution vector (first 5 entries): {solution_k6_float[:5]}")

        # Convert to mpmath
        solution_k6 = np.array([mp.mpc(c) for c in solution_k6_float])

        # Compute A values for K=6
        sqrt2 = mp.sqrt(mp.mpf(2))
        A_k6 = [sqrt2 * mp.sin(L * z / mp.mpf(2)) for z in z_powers_of_2]
        A_conj_k6 = [mp.conj(A) for A in A_k6]

        # Evaluate cubic/quartic sum for K=6
        result_k6 = evaluate_cubic_quartic_sum(K_test, n_test, N, L, z_powers_of_2, solution_k6, A_k6, A_conj_k6)

        print(f"\nK=6 cubic/quartic sum with hsum solution:")
        print(f"  {mp.nstr(result_k6, 20)}")
        print(f"  |result| = {mp.nstr(mp.fabs(result_k6), 15)}")

        if mp.fabs(result_k6) < mp.mpf(1e-20):
            print(f"  ✓ K=6 WORKS: result ≈ 0 (K-induction identity holds on this family!)")
        else:
            print(f"  ✗ K=6 FAILS: result ≠ 0 (K-induction blocked on this family)")

    print(f"\n{'='*120}")
    print("CONCLUSION")
    print(f"{'='*120}")
    if k5_works:
        print(f"Method verified on K=5. K=6 result is trustworthy.")
        if solution_k6_float is not None:
            if mp.fabs(result_k6) < mp.mpf(1e-20):
                print(f"\nK-INDUCTION WORKS on Powers-of-2 family.")
            else:
                print(f"\nK-INDUCTION BLOCKED on Powers-of-2 family (polynomial ≠ 0).")
    else:
        print(f"Method FAILED on K=5 (control test). K=6 result is UNRELIABLE.")

    print(f"{'='*120}")

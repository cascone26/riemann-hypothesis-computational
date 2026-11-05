"""
LAYER 4 COMPLETE TEST: Exact symbolic evaluation of K-induction identity

1. Build hsum matrix with exact entries (num_basis=2 for speed)
2. Solve for exact nullspace solution (via sympy.nullspace())
3. Build cubic/quartic polynomial expressions from CrossK6Key.lean
4. Evaluate polynomial with extracted coefficients
5. Check if result = 0 (symbolically)

Test both K=5 (control) and K=6 (test) on same family.
"""

from sympy import (
    symbols, sqrt, sin, cos, I, Rational, Matrix, legendre, simplify,
    summation, expand, chebyshevu, conjugate, pi as sym_pi
)
import sympy as sp

def vPlusLC_sym(n, L, z):
    """vPlus_{LC}(n) = sqrt(2)*sin(L*z/2) * P_n(z)"""
    sqrt2 = sqrt(2)
    sin_term = sin(L * z / 2)
    leg_n = legendre(n, z)
    return sqrt2 * sin_term * leg_n

def vMinusLC_sym(n, L, z):
    """vMinus_{LC}(n) = sqrt(2)*sin(L*z/2) * U_n(z)"""
    sqrt2 = sqrt(2)
    sin_term = sin(L * z / 2)
    cheb_n = chebyshevu(n, z)
    return sqrt2 * sin_term * cheb_n

def build_hsum_matrix_symbolic(K, z_values, L, num_basis=2):
    """Build hsum matrix with exact symbolic entries."""
    num_rows = num_basis + 1
    num_cols = 4 * K
    M_list = []

    for n in range(num_basis + 1):
        row = []
        for i in range(K):
            z_i = z_values[i]
            z_i_conj = conjugate(z_i)

            v_plus_i = vPlusLC_sym(n, L, z_i)
            v_minus_i = vMinusLC_sym(n, L, z_i)
            v_plus_conj_i = vPlusLC_sym(n, L, z_i_conj)
            v_minus_conj_i = vMinusLC_sym(n, L, z_i_conj)

            row.extend([v_plus_i, v_minus_i, v_plus_conj_i, v_minus_conj_i])

        M_list.append(row)

    return Matrix(M_list)

def find_nullspace_solution(M):
    """Find exact nullspace and return first basis vector."""
    NS = M.nullspace()
    if len(NS) == 0:
        return None
    return NS[0]

def evaluate_cubic_quartic_sum(K, n_var, N, L, z_values, coeffs):
    """
    Evaluate the cubic/quartic sum from CrossK6Key.lean.

    cubic_i = (p + q)*A*z*(conj(z))^2 + (r+s)*conj(A)*z^2*conj(z) + ...
    quartic_i = (z + xi)*(z - xi)*(conj(z) + xi)*(conj(z) - xi)

    sum = Σ cubic_i * (∏_{j≠i} quartic_j)
    """
    sqrt2 = sqrt(2)
    xi = 2 * sym_pi * (n_var - N) / L

    # Extract coefficients
    p_vals = [coeffs[4*i+0] for i in range(K)]
    q_vals = [coeffs[4*i+1] for i in range(K)]
    r_vals = [coeffs[4*i+2] for i in range(K)]
    s_vals = [coeffs[4*i+3] for i in range(K)]

    # Compute A values
    A_vals = [sqrt2 * sin(L * z / 2) for z in z_values]
    A_conj_vals = [conjugate(A) for A in A_vals]

    # Compute quartics
    quartics = []
    for i in range(K):
        z_i = z_values[i]
        z_i_conj = conjugate(z_i)
        q_i = (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)
        quartics.append(q_i)

    # Compute cubics
    cubics = []
    for i in range(K):
        z_i = z_values[i]
        z_i_conj = conjugate(z_i)
        A_i = A_vals[i]
        A_i_conj = A_conj_vals[i]
        p_i, q_i, r_i, s_i = p_vals[i], q_vals[i], r_vals[i], s_vals[i]

        xi2 = xi ** 2
        xi3 = xi ** 3

        # cubic_i construction
        c_i = (
            (p_i + q_i) * A_i * z_i * (z_i_conj ** 2)
            + (r_i + s_i) * A_i_conj * (z_i ** 2) * z_i_conj
            + (q_i - p_i) * A_i * (z_i_conj ** 2) * xi
            + (s_i - r_i) * A_i_conj * (z_i ** 2) * xi
            - (p_i + q_i) * A_i * z_i * xi2
            - (r_i + s_i) * A_i_conj * z_i_conj * xi2
            + (p_i - q_i) * A_i * xi3
            + (r_i - s_i) * A_i_conj * xi3
        )
        cubics.append(c_i)

    # Sum: cubic_i * (all other quartics)
    total = 0
    for i in range(K):
        product_others = 1
        for j in range(K):
            if i != j:
                product_others *= quartics[j]
        total += cubics[i] * product_others

    return total

if __name__ == "__main__":
    print("=" * 120)
    print("LAYER 4: COMPLETE K-INDUCTION TEST (Exact Symbolic)")
    print("=" * 120)

    # Exact family: [0, 1, 2, 3, 4, 5] + (5/2)j
    z_family = [i + Rational(5, 2) * I for i in range(6)]

    L = 20
    N = 10
    n_test = 15
    num_basis = 2

    print(f"\nFamily: [0, 1, 2, 3, 4, 5] + (5/2)j (exact Rational)")
    print(f"L = {L}, N = {N}, test at n = {n_test}")
    print(f"num_basis = {num_basis}")

    # ========== K=5 CONTROL ==========
    print(f"\n{'='*120}")
    print("K=5 CONTROL (known-working)")
    print(f"{'='*120}")

    K = 5
    z_k5 = z_family[:K]

    print(f"\nBuilding hsum matrix for K=5...")
    M_k5 = build_hsum_matrix_symbolic(K, z_k5, L, num_basis=num_basis)
    print(f"Matrix shape: {M_k5.shape}")

    print(f"Finding nullspace solution...")
    x_k5 = find_nullspace_solution(M_k5)

    if x_k5 is None:
        print("ERROR: No nullspace found for K=5")
    else:
        print(f"Solution found (dimension {len(x_k5)})")
        print(f"First 8 entries: {[x_k5[i] for i in range(min(8, len(x_k5)))]}")

        # Evaluate polynomial for K=5
        print(f"\nEvaluating cubic/quartic polynomial for K=5...")
        poly_k5 = evaluate_cubic_quartic_sum(K, n_test, N, L, z_k5, x_k5)

        print(f"Polynomial expression (before simplification):")
        print(f"  (too large to display, length ~ {len(str(poly_k5))} chars)")

        print(f"\nSimplifying...")
        poly_k5_simplified = simplify(poly_k5)

        print(f"Simplified result: {poly_k5_simplified}")

        if poly_k5_simplified == 0:
            print(f"✓ K=5 WORKS: Polynomial = 0 (method validated)")
            k5_works = True
        else:
            print(f"✗ K=5 FAILS: Polynomial ≠ 0 (method bug?)")
            k5_works = False

    # ========== K=6 TEST ==========
    print(f"\n{'='*120}")
    print("K=6 TEST (survivor family)")
    print(f"{'='*120}")

    K = 6
    z_k6 = z_family[:K]

    print(f"\nBuilding hsum matrix for K=6...")
    M_k6 = build_hsum_matrix_symbolic(K, z_k6, L, num_basis=num_basis)
    print(f"Matrix shape: {M_k6.shape}")

    print(f"Finding nullspace solution...")
    x_k6 = find_nullspace_solution(M_k6)

    if x_k6 is None:
        print("ERROR: No nullspace found for K=6")
    else:
        print(f"Solution found (dimension {len(x_k6)})")
        print(f"First 8 entries: {[x_k6[i] for i in range(min(8, len(x_k6)))]}")

        # Evaluate polynomial for K=6
        print(f"\nEvaluating cubic/quartic polynomial for K=6...")
        poly_k6 = evaluate_cubic_quartic_sum(K, n_test, N, L, z_k6, x_k6)

        print(f"Polynomial expression (before simplification):")
        print(f"  (too large to display, length ~ {len(str(poly_k6))} chars)")

        print(f"\nSimplifying (this may take a moment)...")
        poly_k6_simplified = simplify(poly_k6)

        print(f"Simplified result: {poly_k6_simplified}")

        if poly_k6_simplified == 0:
            print(f"✓ K=6 WORKS: Polynomial = 0 (K-induction identity holds!)")
        else:
            print(f"✗ K=6 FAILS: Polynomial ≠ 0 (K-induction blocked)")

    # Final report
    print(f"\n{'='*120}")
    print("FINAL REPORT")
    print(f"{'='*120}")

    if 'k5_works' in locals():
        if k5_works:
            print(f"Method VALIDATED on K=5: symbolic computation reproduces known-working case")
            if x_k6 is not None:
                print(f"K=6 result is TRUSTWORTHY")
                if poly_k6_simplified == 0:
                    print(f"\n*** CONCLUSION: K-INDUCTION WORKS on this family ***")
                else:
                    print(f"\n*** CONCLUSION: K-INDUCTION BLOCKED on this family ***")
        else:
            print(f"Method FAILED on K=5: symbolic computation doesn't reproduce known case")
            print(f"K=6 result is UNRELIABLE")
    else:
        print(f"K=5 test did not complete")

    print(f"{'='*120}")

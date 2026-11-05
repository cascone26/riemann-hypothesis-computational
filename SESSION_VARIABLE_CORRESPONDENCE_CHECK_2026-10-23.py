"""
K=5 CONTROL: Diagnostic version

CRITICAL CHECK: Is the hsum constraint ACTUALLY satisfied?
And does variable order match exactly between hsum and cubic construction?

Test with a GENERIC z-family (not [0,1,2,3,4]+(5/2)j which causes single-block degeneracy)
"""

from sympy import (
    symbols, sqrt, sin, cos, I, Rational, Matrix, legendre, simplify, conjugate,
    chebyshevu, pi as sym_pi, expand, S
)
import sympy as sp

def vPlusLC_sym(n, L, z):
    sqrt2 = sqrt(2)
    sin_term = sin(L * z / 2)
    leg_n = legendre(n, z)
    return sqrt2 * sin_term * leg_n

def vMinusLC_sym(n, L, z):
    sqrt2 = sqrt(2)
    sin_term = sin(L * z / 2)
    cheb_n = chebyshevu(n, z)
    return sqrt2 * sin_term * cheb_n

def build_hsum_matrix(K, z_values, L, n_indices):
    M_list = []
    for n in n_indices:
        row = []
        for i in range(K):
            z_i = z_values[i]
            z_i_conj = conjugate(z_i)
            row.extend([
                vPlusLC_sym(n, L, z_i),
                vMinusLC_sym(n, L, z_i),
                vPlusLC_sym(n, L, z_i_conj),
                vMinusLC_sym(n, L, z_i_conj)
            ])
        M_list.append(row)
    return Matrix(M_list)

def evaluate_hsum_constraint(K, n_eval, L, z_values, coeffs):
    """
    Directly evaluate the hsum constraint at a specific n.

    From CrossK5Key.lean lines 18-27:
      p1*vPlus(z1) + q1*vMinus(z1) + r1*vPlus(conj(z1)) + s1*vMinus(conj(z1))
      + p2*vPlus(z2) + q2*vMinus(z2) + r2*vPlus(conj(z2)) + s2*vMinus(conj(z2))
      + ... = 0
    """
    total = 0
    for i in range(K):
        z_i = z_values[i]
        z_i_conj = conjugate(z_i)

        p_i = coeffs[4*i+0]
        q_i = coeffs[4*i+1]
        r_i = coeffs[4*i+2]
        s_i = coeffs[4*i+3]

        term = (
            p_i * vPlusLC_sym(n_eval, L, z_i)
            + q_i * vMinusLC_sym(n_eval, L, z_i)
            + r_i * vPlusLC_sym(n_eval, L, z_i_conj)
            + s_i * vMinusLC_sym(n_eval, L, z_i_conj)
        )
        total += term

    return total

def evaluate_polynomial(K, n_eval, N, L, z_values, coeffs):
    sqrt2 = sqrt(2)
    tau = 2 * sym_pi * (n_eval - N) / L
    tau2 = tau ** 2
    tau3 = tau ** 3

    p_vals = [coeffs[4*i+0] for i in range(K)]
    q_vals = [coeffs[4*i+1] for i in range(K)]
    r_vals = [coeffs[4*i+2] for i in range(K)]
    s_vals = [coeffs[4*i+3] for i in range(K)]

    A_vals = [sqrt2 * sin(L * z / 2) for z in z_values]
    A_conj_vals = [conjugate(A) for A in A_vals]

    quartics = []
    for i in range(K):
        z_i = z_values[i]
        z_i_conj = conjugate(z_i)
        q_i = (z_i + tau) * (z_i - tau) * (z_i_conj + tau) * (z_i_conj - tau)
        quartics.append(q_i)

    cubics = []
    for i in range(K):
        z_i = z_values[i]
        z_i_conj = conjugate(z_i)
        A_i = A_vals[i]
        A_i_conj = A_conj_vals[i]
        p_i, q_i_coeff, r_i, s_i = p_vals[i], q_vals[i], r_vals[i], s_vals[i]

        cubic_i = (
            (p_i + q_i_coeff) * A_i * z_i * (z_i_conj ** 2)
            + (r_i + s_i) * A_i_conj * (z_i ** 2) * z_i_conj
            + (q_i_coeff - p_i) * A_i * (z_i_conj ** 2) * tau
            + (s_i - r_i) * A_i_conj * (z_i ** 2) * tau
            - (p_i + q_i_coeff) * A_i * z_i * tau2
            - (r_i + s_i) * A_i_conj * z_i_conj * tau2
            + (p_i - q_i_coeff) * A_i * tau3
            + (r_i - s_i) * A_i_conj * tau3
        )
        cubics.append(cubic_i)

    total = 0
    for i in range(K):
        product_others = 1
        for j in range(K):
            if i != j:
                product_others *= quartics[j]
        total += cubics[i] * product_others

    return total


def test_family(family_name, z_family, L, N, n_constraint_indices, n_eval):
    print(f"\n{'='*130}")
    print(f"Testing family: {family_name}")
    print(f"{'='*130}")

    K = len(z_family)
    print(f"z values: {z_family}")
    print(f"L={L}, N={N}")
    print(f"Constraint matrix: n ∈ {n_constraint_indices}")
    print(f"Polynomial eval: n = {n_eval}")

    # Build and solve hsum
    M = build_hsum_matrix(K, z_family, L, n_constraint_indices)
    NS = M.nullspace()

    if len(NS) == 0:
        print("  ERROR: No nullspace!")
        return False

    x_sol = NS[0]
    print(f"\n  Nullspace dimension: {len(NS)}")
    print(f"  Solution (first 8): {[x_sol[i] for i in range(min(8, len(x_sol)))]}")

    # DIAGNOSTIC 1: Verify hsum constraint
    hsum_val = evaluate_hsum_constraint(K, n_eval, L, z_family, x_sol)
    hsum_simp = simplify(hsum_val)
    print(f"\n  [DIAGNOSTIC 1] hsum at n={n_eval}: {hsum_simp}")
    if hsum_simp != 0:
        print(f"    WARNING: hsum NOT zero! Constraint not satisfied!")
    else:
        print(f"    ✓ hsum = 0 (constraint satisfied)")

    # DIAGNOSTIC 2: Evaluate polynomial
    poly = evaluate_polynomial(K, n_eval, N, L, z_family, x_sol)
    poly_simp = simplify(poly)
    print(f"\n  [DIAGNOSTIC 2] Polynomial at n={n_eval}: {str(poly_simp)[:100]}{'...' if len(str(poly_simp)) > 100 else ''}")

    if poly_simp == 0:
        print(f"    ✓✓✓ PASS: Polynomial = 0")
        return True
    else:
        print(f"    ✗✗✗ FAIL: Polynomial ≠ 0")
        return False


def main():
    print("="*130)
    print("K=5 CONTROL: DIAGNOSTIC MODE")
    print("="*130)

    L = 20
    N = 10
    n_constraint = [0, 1, 2]
    n_eval = 0

    # Test 1: Original family (causes single-block degeneracy)
    print("\nTest 1: [0,1,2,3,4] + (5/2)j (original, may be degenerate)")
    family1 = [Rational(i) + Rational(5,2)*I for i in range(5)]
    result1 = test_family("family1", family1, L, N, n_constraint, n_eval)

    # Test 2: More generic family
    print("\nTest 2: [1, 1.5, 2, 2.5, 3] + 0.7j (generic)")
    family2 = [Rational(1) + Rational(7,10)*I,
               Rational(3,2) + Rational(7,10)*I,
               Rational(2) + Rational(7,10)*I,
               Rational(5,2) + Rational(7,10)*I,
               Rational(3) + Rational(7,10)*I]
    result2 = test_family("family2", family2, L, N, n_constraint, n_eval)

    # Summary
    print(f"\n{'='*130}")
    print("SUMMARY")
    print(f"{'='*130}")
    print(f"Family 1 (original): {'PASS' if result1 else 'FAIL'}")
    print(f"Family 2 (generic):  {'PASS' if result2 else 'FAIL'}")

    if result1 and result2:
        print("\n✓ Both families pass - variable mapping is CORRECT")
        return True
    elif result2 and not result1:
        print("\n⚠ Family 2 passes but Family 1 fails - Family 1 may be degenerate edge case")
        return True
    else:
        print("\n✗ Both fail - systematic bug in polynomial evaluation or variable mapping")
        return False


if __name__ == "__main__":
    result = main()
    exit(0 if result else 1)

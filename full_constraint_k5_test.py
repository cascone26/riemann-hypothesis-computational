"""
COMPLETE K=5 CONTROL TEST: Full hsum matrix with N=3 (7 equations × 20 unknowns)

Reduced from N=10 (21 equations) to N=3 (7 equations) for computational tractability.
Still over-constrained (7 > 20 is impossible, meant to say 7 eqs but 20 unknowns → underdetermined).

Actually, let me reconsider: 7 equations, 20 unknowns → underdetermined, has non-trivial nullspace.
10 equations would be more interesting: 10 equations, 20 unknowns.

Let me use N=4, which gives 2*4+1 = 9 equations × 20 unknowns.

Key question: if we increase constraint count by using more n values, does nullspace shrink?
And does any nontrivial solution, if it exists, still give polynomial=0?
"""

from sympy import (
    sqrt, sin, I, Rational, Matrix, legendre, simplify, conjugate,
    chebyshevu, pi as sym_pi, symbols, expand, nsimplify
)
import time

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

def build_complete_hsum_matrix(K, z_values, L, N):
    """Build hsum matrix for all n ∈ {0, ..., 2N}."""
    num_rows = 2 * N + 1
    num_cols = 4 * K
    M_list = []

    for n in range(num_rows):
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

def evaluate_polynomial(K, n_eval, N, L, z_values, coeffs):
    """Evaluate cubic/quartic polynomial."""
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


def test_n_value(N, K, z_family, L):
    """Test with specific N value."""
    num_eqs = 2 * N + 1

    print(f"\n{'='*130}")
    print(f"Testing with N={N}: {num_eqs} equations × 20 unknowns")
    print(f"{'='*130}")

    # Build matrix
    print(f"Building {num_eqs}×20 matrix...")
    t_start = time.time()
    try:
        M = build_complete_hsum_matrix(K, z_family, L, N)
        t_build = time.time() - t_start
        print(f"  Built in {t_build:.2f}s")
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

    # Compute nullspace
    print(f"Computing nullspace...")
    t_start = time.time()
    try:
        NS = M.nullspace()
        t_ns = time.time() - t_start
        print(f"  Nullspace computed in {t_ns:.2f}s")
        print(f"  Dimension: {len(NS)}")
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

    if len(NS) == 0:
        print(f"  Result: TRIVIAL nullspace (all zeros forced)")
        return True
    else:
        print(f"  Result: Nontrivial nullspace (dimension {len(NS)})")
        x_sol = NS[0]

        # Verify
        residual = M @ x_sol
        if all(simplify(r) == 0 for r in residual):
            print(f"  Constraint verified: M @ x = 0 ✓")

            # Evaluate polynomial
            print(f"  Evaluating polynomial...")
            poly = evaluate_polynomial(K, 0, N, L, z_family, x_sol)
            poly_simp = simplify(poly)

            if poly_simp == 0:
                print(f"  Polynomial: 0 ✓")
                return True
            else:
                print(f"  Polynomial: {str(poly_simp)[:80]}... ✗")
                return False
        else:
            print(f"  ERROR: Constraint not satisfied")
            return None


def main():
    print("="*130)
    print("COMPLETE K=5 CONTROL: Graduated constraint scaling")
    print("="*130)

    K = 5
    z_family = [Rational(i) + Rational(5, 2) * I for i in range(K)]
    L = 20

    print(f"Family: {z_family}")
    print(f"Testing with increasing constraint counts...")

    results = {}
    for N in [3, 4, 5]:
        result = test_n_value(N, K, z_family, L)
        results[N] = result

    # Summary
    print(f"\n{'='*130}")
    print("SUMMARY")
    print(f"{'='*130}")
    for N, result in results.items():
        num_eqs = 2*N + 1
        status = "PASS" if result is True else ("FAIL" if result is False else "ERROR")
        print(f"  N={N} ({num_eqs:2d} equations): {status}")

    print(f"\nKey finding:")
    if all(r is True for r in results.values()):
        print(f"  All tests PASS: K=5 control verified at multiple constraint levels")
        return True
    elif any(r is None for r in results.values()):
        print(f"  Some tests errored: computation may be too expensive")
        return None
    else:
        print(f"  Tests show mixed results")
        return False


if __name__ == "__main__":
    result = main()
    exit(0 if result is True else 1)

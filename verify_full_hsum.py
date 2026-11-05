"""
CRITICAL CHECK: Does the partial hsum solution (constrained for n ∈ {0,1,2})
also satisfy the FULL hsum constraint for ALL n ∈ {0,...,2N}?

This is the key issue: the theorem requires hsum to hold for ALL n.
But we're only constraining it for 3 values of n.
"""

from sympy import (
    sqrt, sin, I, Rational, Matrix, legendre, simplify, conjugate,
    chebyshevu, pi as sym_pi, symbols, N as sym_N
)

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

def evaluate_hsum_at_n(K, n_eval, L, z_values, coeffs):
    """Evaluate full hsum constraint at a specific n."""
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

def main():
    print("="*130)
    print("CRITICAL CHECK: Full hsum constraint verification")
    print("="*130)

    K = 5
    z_family = [Rational(i) + Rational(5,2)*I for i in range(K)]
    L = 20
    N = 10

    # Build partial hsum matrix (n ∈ {0, 1, 2})
    print(f"\n[1] Building partial hsum matrix (n ∈ {{0, 1, 2}})...")
    M_partial = build_hsum_matrix(K, z_family, L, [0, 1, 2])
    NS = M_partial.nullspace()
    x_sol = NS[0]
    print(f"    Nullspace dimension: {len(NS)}")
    print(f"    Solution (first 8): {[x_sol[i] for i in range(min(8, len(x_sol)))]}")

    # Now evaluate full hsum constraint at ALL n
    print(f"\n[2] Evaluating FULL hsum constraint at each n ∈ {{0, 1, ..., 2N}}...")
    print(f"    (This shows whether partial solution satisfies full constraint)")

    violation_count = 0
    for n_test in range(2*N + 1):
        hsum_val = evaluate_hsum_at_n(K, n_test, L, z_family, x_sol)
        hsum_simp = simplify(hsum_val)

        # Only print a sample (too verbose otherwise)
        if n_test <= 5 or n_test >= 2*N-1 or hsum_simp != 0:
            status = "✓ SATISFIED" if hsum_simp == 0 else "✗ VIOLATED"
            expr_str = str(hsum_simp)[:60] if hsum_simp != 0 else "0"
            print(f"    n={n_test:2d}: {status} -- {expr_str}{'...' if len(str(hsum_simp)) > 60 else ''}")

        if hsum_simp != 0:
            violation_count += 1

    print(f"\n[3] Summary:")
    print(f"    Constraint satisfied at {2*N+1 - violation_count} / {2*N+1} values of n")
    print(f"    Violations at {violation_count} points")

    print(f"\n{'='*130}")
    if violation_count == 0:
        print("✓ FULL CONSTRAINT SATISFIED")
        print("The partial solution actually satisfies hsum for ALL n!")
        print("Variable mapping is likely correct. Polynomial evaluation bug remains.")
    else:
        print(f"✗ CONSTRAINT VIOLATED at {violation_count} points")
        print("The solution only satisfies hsum for n ∈ {0,1,2}, NOT for all n.")
        print("This explains the polynomial failure: theorem requires hsum for ALL n!")
        print("To properly test, we need to solve the FULL 21×20 system.")
    print(f"{'='*130}")


if __name__ == "__main__":
    main()

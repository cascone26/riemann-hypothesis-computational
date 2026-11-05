"""
FINAL ATTEMPT: Use irregular (non-symmetric) family to force multi-block coupling.

Test family: DELIBERATELY chosen to have NO arithmetic progression, NO symmetry in Re or Im parts.
[3/10+21*I/10, 17/10+23*I/10, 9/10+26*I/10, 12/5+11*I/5, 11/10+14*I/5]
= [0.3+2.1j, 1.7+2.3j, 0.9+2.6j, 2.4+2.2j, 1.1+2.8j] (exact Rational representation)

Expected: nullspace solution has MULTIPLE nonzero blocks (not just one), enabling meaningful cross-term evaluation.
"""

from sympy import (
    symbols, sqrt, sin, cos, I, Rational, Matrix, legendre, simplify, conjugate,
    chebyshevu, pi as sym_pi
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

def build_hsum_k5(z_vals, L, num_basis=2):
    """Build hsum matrix for K=5."""
    num_rows = num_basis + 1
    M_list = []
    for n in range(num_basis + 1):
        row = []
        for i in range(5):
            z_i = z_vals[i]
            z_i_conj = conjugate(z_i)
            row.extend([
                vPlusLC_sym(n, L, z_i),
                vMinusLC_sym(n, L, z_i),
                vPlusLC_sym(n, L, z_i_conj),
                vMinusLC_sym(n, L, z_i_conj)
            ])
        M_list.append(row)
    return Matrix(M_list)

def evaluate_k5_polynomial_literal(K, n_val, N, L, z_vals, p_vals, q_vals, r_vals, s_vals):
    """Evaluate K=5 polynomial."""
    sqrt2 = sqrt(2)
    xi = 2 * sym_pi * (n_val - N) / L

    A_vals = [sqrt2 * sin(L * z / 2) for z in z_vals]
    A_conj_vals = [conjugate(A) for A in A_vals]

    cubics = []
    for i in range(K):
        z_i = z_vals[i]
        z_i_conj = conjugate(z_i)
        A_i = A_vals[i]
        A_i_conj = A_conj_vals[i]
        p_i, q_i, r_i, s_i = p_vals[i], q_vals[i], r_vals[i], s_vals[i]

        cubic_i = (
            (p_i + q_i) * A_i * (z_i * (z_i_conj ** 2))
            + (r_i + s_i) * A_i_conj * ((z_i ** 2) * z_i_conj)
            + ((q_i - p_i) * A_i * (z_i_conj ** 2) + (s_i - r_i) * A_i_conj * (z_i ** 2)) * xi
            + (-(p_i + q_i) * A_i * z_i - (r_i + s_i) * A_i_conj * z_i_conj) * (xi ** 2)
            + ((p_i - q_i) * A_i + (r_i - s_i) * A_i_conj) * (xi ** 3)
        )
        cubics.append(cubic_i)

    quartics = []
    for i in range(K):
        z_i = z_vals[i]
        z_i_conj = conjugate(z_i)
        quartic_i = (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)
        quartics.append(quartic_i)

    total = (
        cubics[0] * quartics[1] * quartics[2] * quartics[3] * quartics[4]
        + cubics[1] * quartics[0] * quartics[2] * quartics[3] * quartics[4]
        + cubics[2] * quartics[0] * quartics[1] * quartics[3] * quartics[4]
        + cubics[3] * quartics[0] * quartics[1] * quartics[2] * quartics[4]
        + cubics[4] * quartics[0] * quartics[1] * quartics[2] * quartics[3]
    )

    return total

if __name__ == "__main__":
    print("=" * 120)
    print("K=5 FINAL ATTEMPT: IRREGULAR (NON-SYMMETRIC) FAMILY")
    print("=" * 120)

    # Irregular family (no symmetry, no arithmetic progression)
    z_family = [
        Rational(3, 10) + Rational(21, 10) * I,
        Rational(17, 10) + Rational(23, 10) * I,
        Rational(9, 10) + Rational(26, 10) * I,
        Rational(12, 5) + Rational(11, 5) * I,
        Rational(11, 10) + Rational(14, 5) * I
    ]

    K = 5
    L = 20
    N = 10
    n_test = 15
    num_basis = 2

    print(f"\nFamily (irregular, no symmetry):")
    for i, z in enumerate(z_family):
        print(f"  z[{i}] = {z}")

    print(f"\nK = {K}, L = {L}, N = {N}, n = {n_test}, num_basis = {num_basis}")

    # Build hsum matrix
    print(f"\nBuilding hsum matrix for irregular family...")
    M = build_hsum_k5(z_family, L, num_basis=num_basis)
    print(f"Matrix shape: {M.shape}")

    # Find nullspace
    print(f"Finding nullspace...")
    NS = M.nullspace()
    if len(NS) == 0:
        print("ERROR: No nullspace found")
    else:
        x = NS[0]
        print(f"Solution found (dimension {len(x)})")

        # Check for multi-block activity
        print(f"\nAnalyzing solution structure:")
        nonzero_blocks = []
        for block_idx in range(K):
            p, q, r, s = x[4*block_idx], x[4*block_idx+1], x[4*block_idx+2], x[4*block_idx+3]
            is_nonzero = (p != 0) or (q != 0) or (r != 0) or (s != 0)
            if is_nonzero:
                nonzero_blocks.append(block_idx)
                print(f"  Block {block_idx}: p={p}, q={q}, r={r}, s={s}")
            else:
                print(f"  Block {block_idx}: [all zero]")

        if len(nonzero_blocks) == 1:
            print(f"\n✗ DEGENERATE: Only block {nonzero_blocks[0]} is active (single-block coupling)")
            print(f"  This family also decouples — cannot test cross-block interaction")
        elif len(nonzero_blocks) >= 2:
            print(f"\n✓ MULTI-BLOCK: {len(nonzero_blocks)} blocks active — meaningful cross-term test possible")

            # Extract coefficients
            p_vals = [x[4*i] for i in range(K)]
            q_vals = [x[4*i+1] for i in range(K)]
            r_vals = [x[4*i+2] for i in range(K)]
            s_vals = [x[4*i+3] for i in range(K)]

            # Evaluate polynomial
            print(f"\nEvaluating K=5 polynomial with multi-block solution...")
            poly = evaluate_k5_polynomial_literal(K, n_test, N, L, z_family, p_vals, q_vals, r_vals, s_vals)

            print(f"Simplifying (this may take a moment)...")
            poly_simplified = simplify(poly)

            print(f"\nSimplified result: {poly_simplified}")

            if poly_simplified == 0:
                print(f"\n✓✓ K=5 CONTROL PASSES: Polynomial = 0 (MULTI-BLOCK SOLUTION)")
                print(f"METHOD VALIDATED — proceeding to K=6 test with same family")
                k5_passes = True
            else:
                print(f"\n✗ K=5 STILL FAILS: Polynomial ≠ 0 (even with multi-block solution)")
                k5_passes = False
        else:
            print(f"\nERROR: No nonzero blocks in solution")
            k5_passes = False

    print(f"\n{'='*120}")
    if 'k5_passes' in locals() and k5_passes:
        print("READY FOR K=6 TEST")
    else:
        print("UNABLE TO PROCEED (K=5 control still fails)")
    print(f"{'='*120}")

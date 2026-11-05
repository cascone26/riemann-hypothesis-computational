"""
LITERAL TRANSCRIPTION of K=5 cubic from CrossK5Key.lean lines 34-38.

cubic1 := (((p1 + q1) * A1 * (z1 * (starRingEnd ℂ z1) ^ 2)
    + (r1 + s1) * (starRingEnd ℂ A1) * (z1 ^ 2 * (starRingEnd ℂ z1))
    + ((q1 - p1) * A1 * (starRingEnd ℂ z1) ^ 2 + (s1 - r1) * (starRingEnd ℂ A1) * z1 ^ 2) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ)
    + (-(p1 + q1) * A1 * z1 - (r1 + s1) * (starRingEnd ℂ A1) * (starRingEnd ℂ z1)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 2
    + ((p1 - q1) * A1 + (r1 - s1) * (starRingEnd ℂ A1)) * (2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℂ) ^ 3))

Transcribed to sympy:
- starRingEnd ℂ z = conjugate(z)
- A = sqrt(2) * sin(L*z/2)
- xi = 2*pi*(n - N)/L
"""

from sympy import (
    symbols, sqrt, sin, cos, I, Rational, Matrix, legendre, simplify, conjugate,
    chebyshevu, pi as sym_pi, expand
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
    """
    Evaluate K=5 polynomial EXACTLY as written in CrossK5Key.lean lines 34-38, 59-63, 64-68.

    cubic_i = (p_i + q_i) * A_i * z_i * conj(z_i)^2
            + (r_i + s_i) * conj(A_i) * z_i^2 * conj(z_i)
            + ((q_i - p_i) * A_i * conj(z_i)^2 + (s_i - r_i) * conj(A_i) * z_i^2) * xi
            - (p_i + q_i) * A_i * z_i * xi^2
            - (r_i + s_i) * conj(A_i) * conj(z_i) * xi^2
            + (p_i - q_i) * A_i * xi^3
            + (r_i - s_i) * conj(A_i) * xi^3

    quartic_i = (z_i + xi) * (z_i - xi) * (conj(z_i) + xi) * (conj(z_i) - xi)

    sum = cubic_1 * Q_2 * Q_3 * Q_4 * Q_5
        + cubic_2 * Q_1 * Q_3 * Q_4 * Q_5
        + cubic_3 * Q_1 * Q_2 * Q_4 * Q_5
        + cubic_4 * Q_1 * Q_2 * Q_3 * Q_5
        + cubic_5 * Q_1 * Q_2 * Q_3 * Q_4
    """
    sqrt2 = sqrt(2)
    xi = 2 * sym_pi * (n_val - N) / L

    # Compute A_i values
    A_vals = [sqrt2 * sin(L * z / 2) for z in z_vals]
    A_conj_vals = [conjugate(A) for A in A_vals]

    # Compute cubics (LITERAL TRANSCRIPTION)
    cubics = []
    for i in range(K):
        z_i = z_vals[i]
        z_i_conj = conjugate(z_i)
        A_i = A_vals[i]
        A_i_conj = A_conj_vals[i]
        p_i, q_i, r_i, s_i = p_vals[i], q_vals[i], r_vals[i], s_vals[i]

        # Line 34-38 from CrossK5Key.lean:
        cubic_i = (
            (p_i + q_i) * A_i * (z_i * (z_i_conj ** 2))
            + (r_i + s_i) * A_i_conj * ((z_i ** 2) * z_i_conj)
            + ((q_i - p_i) * A_i * (z_i_conj ** 2) + (s_i - r_i) * A_i_conj * (z_i ** 2)) * xi
            + (-(p_i + q_i) * A_i * z_i - (r_i + s_i) * A_i_conj * z_i_conj) * (xi ** 2)
            + ((p_i - q_i) * A_i + (r_i - s_i) * A_i_conj) * (xi ** 3)
        )
        cubics.append(cubic_i)

    # Compute quartics (LINE 59-63 from CrossK5Key.lean)
    quartics = []
    for i in range(K):
        z_i = z_vals[i]
        z_i_conj = conjugate(z_i)
        quartic_i = (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)
        quartics.append(quartic_i)

    # Sum (LINES 64-68 from CrossK5Key.lean)
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
    print("K=5 CONTROL: LITERAL TRANSCRIPTION from CrossK5Key.lean")
    print("=" * 120)

    # Exact family: [0, 1, 2, 3, 4] + (5/2)j
    z_family = [i + Rational(5, 2) * I for i in range(5)]

    K = 5
    L = 20
    N = 10
    n_test = 15
    num_basis = 2

    print(f"\nFamily: [0, 1, 2, 3, 4] + (5/2)j")
    print(f"K = {K}, L = {L}, N = {N}, n = {n_test}")
    print(f"num_basis = {num_basis}")

    # Build hsum matrix
    print(f"\nBuilding hsum matrix...")
    M = build_hsum_k5(z_family, L, num_basis=num_basis)
    print(f"Matrix shape: {M.shape}")

    # Find nullspace
    print(f"Finding nullspace...")
    NS = M.nullspace()
    if len(NS) == 0:
        print("ERROR: No nullspace")
    else:
        x = NS[0]
        print(f"Solution found (dimension {len(x)})")
        print(f"First 8 entries: {[x[i] for i in range(min(8, len(x)))]}")

        # Extract coefficients for each z_i
        p_vals = [x[4*i+0] for i in range(K)]
        q_vals = [x[4*i+1] for i in range(K)]
        r_vals = [x[4*i+2] for i in range(K)]
        s_vals = [x[4*i+3] for i in range(K)]

        print(f"\nCoefficients:")
        for i in range(K):
            print(f"  z[{i}]: p={p_vals[i]}, q={q_vals[i]}, r={r_vals[i]}, s={s_vals[i]}")

        # Evaluate polynomial
        print(f"\nEvaluating K=5 polynomial (LITERAL TRANSCRIPTION)...")
        poly = evaluate_k5_polynomial_literal(K, n_test, N, L, z_family, p_vals, q_vals, r_vals, s_vals)

        print(f"Polynomial expression (before simplification):")
        print(f"  (length ~ {len(str(poly))} chars)")

        print(f"\nSimplifying...")
        poly_simplified = simplify(poly)

        print(f"Simplified result: {poly_simplified}")

        if poly_simplified == 0:
            print(f"\n✓ K=5 WORKS: Polynomial = 0 (LITERAL TRANSCRIPTION VERIFIED)")
        else:
            print(f"\n✗ K=5 FAILS: Polynomial ≠ 0")
            print(f"Result magnitude: {simplify(sp.Abs(poly_simplified))}")

    print(f"\n{'='*120}")

#!/usr/bin/env python3
"""
Characterize K-dependent rank behavior of vPlusLC/vMinusLC vectors.

The K-induction fails at K=6 because the hsum constraint has K-dependent properties.
The hsum constraint is: sum_i [p_i*v+(z_i) + q_i*v-(z_i) + r_i*v+(conj z_i) + s_i*v-(conj z_i)] = 0

This constraint defines a K-dependent solution space. For K=3,4,5, this space has
cancellation properties that force difference terms (p-q) and (r-s) to vanish. For K=6,
they don't. The root cause is the rank structure of the 4K vector family.

Task: compute the rank for K=3,4,5,6 and identify the exact K=6 dependency.
"""

import numpy as np
import sympy as sp
from sympy import symbols, sin, cos, I, re, im, conjugate, simplify, expand, Matrix, zeros

def vPlusLC_symbolic(N, L, z):
    """
    vPlusLC symbolic formula. From the code and literature:
    v+(z) = sin(Lz/2) / sin(L/2) for the basic form, but scaled differently in actual implementation.

    Looking at the Lean code and gen_cross_k_key.py, vPlusLC involves:
    A(z) = sqrt(2) * sin(Lz/2)

    The actual formula for vPlusLC appears to be a vector-valued function indexed by (N,L,z),
    with N+1 components (one per frequency bin n=0..2N).

    We'll use a simplified symbolic representation based on the sinc structure:
    v+(z)[n] ≈ sin(Lz/2) / sin(L/2) * kernel(z, n)

    But for rank analysis, we can work with the abstract structure: vPlusLC is a function
    that takes (N, L, z) -> vector in C^(2N+1).
    """
    # Simplified: represent vPlusLC as a symbolic function of z, parametrized by (N, L)
    # The key is that v+(z), v-(z), v+(conj z), v-(conj z) form a 4-vector family per z
    # For rank analysis, we focus on the functional dependence structure

    # Using the sin(Lz/2) factor as the core:
    A = sin(L * z / 2)
    # Normalize by sin(L/2) for consistency
    return A / sin(L / 2) if sin(L / 2) != 0 else A

def vMinusLC_symbolic(N, L, z):
    """
    vMinusLC symbolic formula. Similar to vPlusLC but with a different phase/structure.
    From the context, vMinusLC appears to be conjugate-related or phase-shifted.
    """
    A = sin(L * z / 2)
    return A.conjugate() / sin(L / 2) if sin(L / 2) != 0 else A.conjugate()

def build_hsum_constraint_symbolic(K, N, L, use_conjugate=True):
    """
    Build the hsum constraint symbolically for K parameters z1..zK.

    hsum = sum_i [p_i*v+(z_i) + q_i*v-(z_i) + r_i*v+(conj z_i) + s_i*v-(conj z_i)] = 0

    This is a linear constraint on the 4K coefficients (p1,q1,r1,s1,...,pK,qK,rK,sK).
    The constraint is indexed by n (frequency bin), so we have 2N+1 constraints total.

    For rank analysis, we extract the structure: the constraint matrix has rows indexed by n,
    columns indexed by the 4K coefficients, and entries are evaluations of the 4K vectors.
    """

    # Symbolic parameters
    zs = symbols(f'z1:{K+1}')
    ps = symbols(f'p1:{K+1}')
    qs = symbols(f'q1:{K+1}')
    rs = symbols(f'r1:{K+1}')
    ss = symbols(f's1:{K+1}')

    # For K-dependent rank analysis, we examine the structure of the coefficient vectors.
    # At a given frequency n, the hsum constraint involves:
    # [v+(z_1)[n], v-(z_1)[n], v+(conj z_1)[n], v-(conj z_1)[n], ..., v+(z_K)[n], ..., v-(conj z_K)[n]]
    # multiplied by the coefficient vector [p_1, q_1, r_1, s_1, ..., p_K, q_K, r_K, s_K]

    # For symbolic analysis, we focus on one specific frequency n (say, n=0 or a generic n)
    # and analyze when the 4K vectors are linearly independent.

    # Build the vector family for a single frequency evaluation:
    vectors = []
    for i in range(K):
        # v+(z_i), v-(z_i), v+(conj z_i), v-(conj z_i)
        v_plus = vPlusLC_symbolic(N, L, zs[i])
        v_minus = vMinusLC_symbolic(N, L, zs[i])
        v_plus_conj = vPlusLC_symbolic(N, L, conjugate(zs[i]))
        v_minus_conj = vMinusLC_symbolic(N, L, conjugate(zs[i]))

        vectors.extend([v_plus, v_minus, v_plus_conj, v_minus_conj])

    return vectors, zs, ps, qs, rs, ss

def analyze_rank_numerical(K, N=4, L_val=2*np.pi):
    """
    Numerical rank analysis for K parameters with random z values.

    We'll compute the Gram matrix of the 4K vectors evaluated at random points
    to infer the rank structure.
    """
    np.random.seed(42 + K)  # Reproducible randomness per K

    # Random z values (off critical line)
    zs = np.array([np.random.randn() + 1j*np.random.uniform(0.1, 2, 1) for _ in range(K)]).flatten()

    # Evaluate the vectors at a specific point (say, n=0 or frequency 0)
    # For simplicity, we compute vPlusLC/vMinusLC values directly

    vectors = []
    for z_val in zs:
        # v+(z) ~ sin(Lz/2) / sin(L/2)
        A = np.sin(L_val * z_val / 2) / np.sin(L_val / 2) if np.sin(L_val / 2) != 0 else np.sin(L_val * z_val / 2)
        A_conj = np.conj(A)

        # For a rank analysis, we evaluate the vector values as scalars (e.g., at n=0)
        # v_plus(z)[0] ≈ A, v_minus(z)[0] ≈ A_conj, etc.
        v_plus = A
        v_minus = np.conj(A)
        v_plus_conj = A  # v+(conj z) evaluated at n=0
        v_minus_conj = np.conj(A)  # v-(conj z) evaluated at n=0

        vectors.extend([v_plus, v_minus, v_plus_conj, v_minus_conj])

    # Form a "Gram-like" matrix (not exactly Gram, but captures linear relationships)
    # We compute pairwise inner products to infer rank
    vectors = np.array(vectors, dtype=complex)

    # Compute the Gram matrix: G[i,j] = conj(vectors[i]) * vectors[j]
    G = np.outer(np.conj(vectors), vectors)

    # Compute rank via SVD
    U, s, Vh = np.linalg.svd(G)
    rank = np.linalg.matrix_rank(G)

    # Identify small singular values (near-dependencies)
    singular_values = s
    tol = 1e-10
    small_sv_count = np.sum(singular_values < tol)

    return {
        'K': K,
        'num_vectors': 4*K,
        'rank': rank,
        'singular_values': singular_values[:10],  # Top 10 SV
        'small_sv_count': small_sv_count,
        'condition_number': singular_values[0] / (singular_values[-1] + 1e-16),
        'zs': zs
    }

def analyze_hsum_constraint_structure(K):
    """
    Direct analysis of the hsum constraint structure for different K.

    The key insight: for K=3,4,5, the constraint space has dimension 4K-1 (one constraint).
    But the cancellation properties differ at K=6.

    We'll investigate: what is special about K=3,4,5 that K=6 lacks?
    """
    print(f"\n{'='*70}")
    print(f"K={K} ANALYSIS: hsum constraint structure")
    print(f"{'='*70}")

    # Theoretical expectation:
    # - 4K coefficient variables (p1,q1,r1,s1, ..., pK,qK,rK,sK)
    # - 1 constraint (the hsum equation sums to 0)
    # - Solution space dimension: 4K - 1

    # The failure at K=6 suggests that at K=6, the LINEAR DEPENDENCE STRUCTURE
    # of the 4K vectors changes. Either:
    # (a) The rank of the vector family changes (becomes < 4K)
    # (b) The null space structure changes

    sol_space_dim = 4*K - 1  # Expected solution space dimension

    print(f"Theoretical solution space dimension: {sol_space_dim}")
    print(f"(4*{K} variables - 1 constraint = {sol_space_dim})")

    # Numerical rank analysis
    rank_info = analyze_rank_numerical(K)
    print(f"\nNumerical Gram-matrix rank: {rank_info['rank']}")
    print(f"Expected rank (full): {4*K}")
    print(f"Singular value spectrum (top 10): {rank_info['singular_values']}")
    print(f"Condition number: {rank_info['condition_number']:.2e}")

    # Check for specific near-dependencies
    if rank_info['small_sv_count'] > 0:
        print(f"\n⚠ WARNING: {rank_info['small_sv_count']} small singular values detected")
        print("  This suggests linear dependencies in the vector family.")
        print("  For K=6, this may be where cancellation breaks down.")

    return rank_info

def main():
    print("\n" + "="*70)
    print("VPLUSLC/VMINUSLC K-DEPENDENT RANK ANALYSIS")
    print("="*70)
    print("\nGoal: characterize why K-induction works for K=3,4,5 but fails at K=6")
    print("Method: analyze the rank structure of the 4K vector family per K")

    results = {}
    for K in [3, 4, 5, 6]:
        results[K] = analyze_hsum_constraint_structure(K)

    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70)
    print(f"{'K':<4} {'Vectors':<10} {'Rank':<8} {'Full?':<8} {'SV_tol':<12} {'Cond#':<12}")
    print("-"*70)
    for K in [3, 4, 5, 6]:
        r = results[K]
        full = "YES" if r['rank'] == 4*K else "NO"
        print(f"{K:<4} {r['num_vectors']:<10} {r['rank']:<8} {full:<8} {r['small_sv_count']:<12} {r['condition_number']:<12.2e}")

    print("\n" + "="*70)
    print("HYPOTHESIS")
    print("="*70)
    print("""
For K=3,4,5: The 4K vectors are linearly independent (or very close).
            The hsum constraint has special cancellation properties that force
            the (p-q) and (r-s) difference terms to vanish in the cubic sum.

For K=6: The rank structure changes. Either:
         (a) The vectors become linearly dependent (rank < 4*6 = 24)
         (b) The null space structure changes, losing the cancellation property

The K=6 failure in the cubic identity is a manifestation of this rank change.
""")

    print("\nNEXT STEP: Identify the exact K=6 dependency or rank deficit")
    print("          by computing the null space of the constraint matrix")
    print("          and comparing it to K=3,4,5.")

if __name__ == "__main__":
    main()

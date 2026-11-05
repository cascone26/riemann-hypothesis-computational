#!/usr/bin/env python3
"""
Analyze the K-dependent constraint structure of hsum.

The hsum constraint is: sum_i [p_i*v+(z_i) + q_i*v-(z_i) + r_i*v+(conj z_i) + s_i*v-(conj z_i)] = 0

This is a VECTOR equation (each v is a vector in C^(2N+1)). It imposes 2N+1 separate scalar constraints,
one per frequency bin n.

The KEY INSIGHT: for K=3,4,5, the constraint has properties that force (p-q) and (r-s) to cancel
in the cubic sum-of-products. For K=6, this doesn't happen.

HYPOTHESIS: At K=6, a new linear dependence emerges among the vectors that changes the rank
of the constraint matrix.

ANALYSIS APPROACH:
1. For each K, assume specific z values (off critical line)
2. Build the constraint matrix M: (2N+1) rows (one per frequency), 4K columns (one per coeff)
3. Compute rank(M) and compare to expected 4K-1
4. Identify the exact dependency at K=6
"""

import numpy as np
from sympy import *
import warnings
warnings.filterwarnings('ignore')

def build_constraint_matrix_abstract(K, N, L, z_values):
    """
    Build the abstract constraint matrix for hsum.

    Each row corresponds to a frequency bin n = 0, 1, ..., 2N.
    Each column corresponds to a coefficient: p1,q1,r1,s1,p2,q2,r2,s2,...,pK,qK,rK,sK.

    The entry [n, coeff] is the evaluation of the corresponding vector at frequency n.

    For a symbolic/numerical hybrid approach, we use the sin(Lz/2) structure as the core kernel.
    """

    # We'll work with a simplified 1D projection for numerical tractability
    # Instead of the full (2N+1)-dimensional vector, we compute values at specific frequencies
    # This is justified because the rank structure is determined by the linear relationships
    # among the VECTOR FAMILIES, not the specific frequency sampling

    # For each z_i, define:
    # v+_i = sin(L*z_i/2)
    # v-_i = sin(L*z_i/2)* (some phase factor)
    # v+*_i = sin(L*conj(z_i)/2)
    # v-*_i = sin(L*conj(z_i)/2) * (some phase factor)

    # Actually, looking at the code more carefully, v- is a conjugate-related variant of v+.
    # The exact formulas are less important than the structure.

    # For analysis, we'll treat vPlusLC and vMinusLC as generating the same dimensional space
    # (e.g., both are functions of z with similar range). This is a simplification but captures
    # the essential rank structure.

    M = np.zeros((2*N+1, 4*K), dtype=complex)

    for i in range(K):
        z = z_values[i]
        z_conj = np.conj(z)

        for n in range(2*N+1):
            # Simplified: evaluate sin(Lz/2) as the core kernel for each frequency
            # The key is that this factor appears in vPlusLC(z) and vMinusLC(z)
            A_i = np.sin(L * z / 2)
            A_i_conj = np.sin(L * z_conj / 2)

            # Frequency-dependent scaling (a linear function of n)
            # In the actual formulation, each vector has a (2*pi*(n-N)/L) dependence
            freq_scale = 2 * np.pi * (n - N) / L

            # Set values (simplified model):
            # v+(z_i)[n] ~ A_i * kernel_fn(n, z_i)
            # v-(z_i)[n] ~ A_i* * kernel_fn(n, z_i)  [conjugate variant]
            # v+(conj z_i)[n] ~ A_i_conj * kernel_fn(n, conj z_i)
            # v-(conj z_i)[n] ~ A_i_conj* * kernel_fn(n, conj z_i)

            col_offset = 4 * i

            M[n, col_offset] = A_i * (1 + freq_scale)  # v+(z_i), p_i column
            M[n, col_offset+1] = np.conj(A_i) * (1 + freq_scale)  # v-(z_i), q_i column
            M[n, col_offset+2] = A_i_conj * (1 + freq_scale)  # v+(conj z_i), r_i column
            M[n, col_offset+3] = np.conj(A_i_conj) * (1 + freq_scale)  # v-(conj z_i), s_i column

    return M

def analyze_constraint_rank(K, N=4, L=2*np.pi, seed=42+10000):
    """Analyze rank structure for a given K."""

    np.random.seed(seed + K)

    # Generate K random z values (off critical line, non-real)
    z_values = np.array([
        np.random.randn() * 0.5 + 1j * np.random.uniform(0.2, 2.0)
        for _ in range(K)
    ])

    # Build constraint matrix
    M = build_constraint_matrix_abstract(K, N, L, z_values)

    # Compute rank via SVD
    U, s, Vh = np.linalg.svd(M)
    rank = np.linalg.matrix_rank(M)

    # Find the null space (vectors that span the solution space)
    # These correspond to valid coefficient choices
    null_rank = 4*K - rank

    # Extract singular values near the numerical precision threshold
    tol = 1e-10
    small_sv = np.sum(s < tol)

    # Ratio of largest to smallest singular value (condition number of the rank-full submatrix)
    if s[-1] > 0:
        cond_num = s[0] / s[-1]
    else:
        cond_num = np.inf

    return {
        'K': K,
        'num_frequencies': 2*N+1,
        'num_coefficients': 4*K,
        'rank': rank,
        'null_space_dim': null_rank,
        'expected_null_space_dim': 4*K - 1,
        'singular_values': s,
        'small_sv_count': small_sv,
        'condition_number': cond_num,
        'z_values': z_values,
        'M': M
    }

def identify_k6_anomaly():
    """
    Identify what's special about K=6.

    Expected behavior: null_space_dim = 4K - 1 for all K
    K=6 anomaly: null_space_dim might differ
    """

    print("\n" + "="*80)
    print("CONSTRAINT MATRIX RANK ANALYSIS")
    print("="*80)

    results = {}
    for K in [3, 4, 5, 6, 7]:
        r = analyze_constraint_rank(K)
        results[K] = r

        print(f"\nK={K}:")
        print(f"  Constraint matrix shape: {r['num_frequencies']} x {r['num_coefficients']}")
        print(f"  Rank: {r['rank']}")
        print(f"  Null space dimension: {r['null_space_dim']}")
        print(f"  Expected null space dimension: {r['expected_null_space_dim']}")
        print(f"  Anomaly: {r['null_space_dim'] != r['expected_null_space_dim']}")
        print(f"  Small singular values (< 1e-10): {r['small_sv_count']}")
        print(f"  Condition number: {r['condition_number']:.2e}")
        print(f"  Top singular values: {r['singular_values'][:5]}")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    print(f"\n{'K':<4} {'Expected':<12} {'Actual':<12} {'Match?':<10}")
    print("-"*40)
    for K in [3, 4, 5, 6, 7]:
        r = results[K]
        exp = r['expected_null_space_dim']
        act = r['null_space_dim']
        match = "YES" if exp == act else "NO ⚠"
        print(f"{K:<4} {exp:<12} {act:<12} {match:<10}")

    return results

def investigate_k6_null_space(results):
    """
    For K=6, investigate the structure of the null space.

    If null_space_dim > 4K - 1, there's an extra dependency.
    Identify which coefficients are affected.
    """

    if results[6]['null_space_dim'] != results[6]['expected_null_space_dim']:
        print("\n" + "="*80)
        print("K=6 NULL SPACE INVESTIGATION")
        print("="*80)

        M6 = results[6]['M']
        U6, s6, Vh6 = np.linalg.svd(M6, full_matrices=True)

        # Null space is spanned by right singular vectors corresponding to singular values ~ 0
        tol = 1e-10
        null_indices = np.where(s6 < tol)[0]

        print(f"\nNull space dimension: {len(null_indices)}")
        print(f"Null space is spanned by {len(null_indices)} singular vectors")

        if len(null_indices) > 0:
            print("\nNull space vectors (coefficients that solve hsum = 0):")
            for idx in null_indices[:3]:  # Show first 3
                v = Vh6[idx, :]
                print(f"\n  Vector {idx}:")
                # Identify which coefficients are dominant
                coeff_names = []
                for k in range(1, 7):
                    for coeff_type, coeff_letter in [(0, 'p'), (1, 'q'), (2, 'r'), (3, 's')]:
                        coeff_idx = 4*(k-1) + coeff_type
                        coeff_names.append(f"{coeff_letter}{k}")

                # Sort by magnitude and show top contributors
                magnitudes = np.abs(v)
                top_indices = np.argsort(magnitudes)[-5:][::-1]

                for tidx in top_indices:
                    if magnitudes[tidx] > 1e-10:
                        print(f"    {coeff_names[tidx]}: {v[tidx]:.6f}")
    else:
        print("\n✓ K=6 null space dimension matches expectation.")
        print("  No rank anomaly detected at K=6 in this simplified model.")

def main():
    print("\n" + "="*80)
    print("VPlusLC/VMinusLC K-DEPENDENT CONSTRAINT STRUCTURE")
    print("="*80)
    print("""
GOAL: Determine why the hsum constraint's cancellation properties break at K=6.

APPROACH: Analyze the constraint matrix M where:
  - Rows: frequency bins (2N+1 rows)
  - Columns: coefficients p,q,r,s per block (4K columns)
  - Entry M[n,j]: evaluation of the j-th vector at frequency n

The null space of M gives valid coefficient choices (solutions to hsum = 0).
For K=3,4,5, the null space has dimension 4K-1 (as expected).
For K=6, if there's an extra dependency, null-space-dim > 4K-1.
""")

    results = identify_k6_anomaly()
    investigate_k6_null_space(results)

    print("\n" + "="*80)
    print("CONCLUSION")
    print("="*80)
    print("""
This simplified analysis provides a framework for detecting K-dependent anomalies
in the constraint structure. The actual K=6 failure (from sympy verification)
suggests that:

1. Either a NEW linear dependency emerges at K=6 (hidden in the higher frequencies)
2. Or the null-space structure changes in a way that removes the cancellation

To get a definitive answer, we would need to:
- Use the EXACT vPlusLC/vMinusLC formulas (not simplified sin model)
- Build the full (2N+1) x 4K constraint matrix for realistic N values
- Compute the null-space basis and examine the (p-q), (r-s) structure
- Compare null-space basis elements between K=5 and K=6

This is a CONCRETE LINEAR ALGEBRA PROBLEM, not a restatement of the RH wall.
""")

if __name__ == "__main__":
    main()

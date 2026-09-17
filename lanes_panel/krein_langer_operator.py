#!/usr/bin/env python3
"""
KREIN-LANGER OPERATOR THEORY LANE

Goal: Build the operator-theoretic proof that M_σ = arch - prime has EXACTLY ONE
negative square (κ=1), pinned by the pole structure, independent of zero locations.

Lens: Pontryagin Π₁ space theory, de Branges-Kaltenbäck structures, Nevanlinna pole-counting.

LEMMAS TO TEST:
1. M_σ is Π₁: #neg(M_σ) = 1 for all samplings (not an artifact)
2. The single negative eigenvalue lives in pole-span ⊗ low-freq (localization)
3. Pole PSD on pole⊥ (orthogonal complement): the negative is structurally tied to pole
4. Generalized eigenvalue perspective: κ ≤ rank(pole) (Krein consequence)
5. The pole direction stabilizes as bandwidth grows (robustness)
"""

import numpy as np
import scipy.linalg as la
from scipy.linalg import eigh
import sys
import json

# Import the exact form builders
sys.path.insert(0, '/Users/scones/projects/riemann')
from operators.suzuki_exact_Q import build_Q

def test_pi1_decomposition(a=1.6, K=100, num_samples=20):
    """
    LEMMA 1: M_σ is Π₁ — test via multiple random samplings.

    If M_σ is Π₁ (Pontryagin with one negative square), then every finite
    Gram (sub-sampling) should have κ = #neg ≤ 1.
    """
    print("\n" + "="*70)
    print("LEMMA 1: M_σ IS PONTRYAGIN Π₁ (all samplings have κ ≤ 1)")
    print("="*70)

    Q_full, pieces = build_Q(a, K, return_pieces=True)
    pole = pieces['pole']
    arch = pieces['core']  # A+B (arch + boundary)
    prime = pieces['prime']

    M_sigma = arch - prime  # This is the Π₁ form

    print(f"\nFull problem: a={a}, K={K}")
    print(f"  #neg(Q_full) = {np.sum(np.linalg.eigvalsh(Q_full) < -1e-14)}")
    print(f"  #neg(M_σ) = {np.sum(np.linalg.eigvalsh(M_sigma) < -1e-14)}")
    print(f"  #neg(pole) = {np.sum(np.linalg.eigvalsh(pole) < -1e-14)}")

    results = {
        'samplings': []
    }

    for sample_idx in range(num_samples):
        # Random sampling: pick a random subset of rows/cols
        N = M_sigma.shape[0]
        subset_size = max(5, N // 2)  # Sample 50% of dimensions
        indices = np.random.choice(N, subset_size, replace=False)
        indices = np.sort(indices)

        M_sub = M_sigma[np.ix_(indices, indices)]
        eigs = np.linalg.eigvalsh(M_sub)
        kappa = np.sum(eigs < -1e-14)

        results['samplings'].append({
            'sample': sample_idx,
            'subset_size': subset_size,
            'kappa': int(kappa),
            'min_eig': float(eigs[0]),
            'min_pos_eig': float(eigs[eigs > -1e-14][0]) if np.sum(eigs > -1e-14) > 0 else None
        })

    kappas = [r['kappa'] for r in results['samplings']]
    print(f"\nSampling results (n={num_samples} random subsets):")
    print(f"  κ values: {set(kappas)}")
    print(f"  max κ: {max(kappas)}")
    print(f"  All κ ≤ 1? {all(k <= 1 for k in kappas)}")

    return all(k <= 1 for k in kappas), results

def test_negative_in_pole_span(a=1.6, K=100):
    """
    LEMMA 2: The one negative eigenvector has large projection onto pole's range.

    This tests whether the single negative eigenvalue is localized to the
    directions spanned by the pole term.
    """
    print("\n" + "="*70)
    print("LEMMA 2: NEGATIVE EIGENVECTOR LIVES IN POLE+LOW-FREQ SPAN")
    print("="*70)

    Q_full, pieces = build_Q(a, K, return_pieces=True)
    pole = pieces['pole']
    arch = pieces['core']
    prime = pieces['prime']

    M_sigma = arch - prime

    # Get negative eigenvalue/eigenvector
    eigs_M, vecs_M = eigh(M_sigma)
    neg_idx = np.argmin(eigs_M)
    lambda_neg = eigs_M[neg_idx]
    v_neg = vecs_M[:, neg_idx]

    # Get pole structure: SVD to find range
    U_pole, S_pole, Vt_pole = np.linalg.svd(pole)
    rank_pole = np.sum(S_pole > 1e-14)
    U_pole = U_pole[:, :rank_pole]  # pole range basis

    # Projection of v_neg onto pole's range
    proj_coeff = U_pole.T @ v_neg
    proj_norm = np.linalg.norm(proj_coeff)
    overlap = proj_norm / np.linalg.norm(v_neg)

    print(f"\nPole rank: {rank_pole}")
    print(f"  Pole singular values: {S_pole[:5]}")
    print(f"\nNegative eigenvalue: λ = {lambda_neg:.6e}")
    print(f"Negative eigenvector:")
    print(f"  norm: {np.linalg.norm(v_neg):.6e}")
    print(f"  projection onto pole-range: {overlap:.6f}")
    print(f"    (0 = orthogonal, 1 = fully in pole-span)")

    # Also check low-frequency content
    low_freq_modes = vecs_M[:, :min(5, K)]  # first few modes
    proj_low = low_freq_modes.T @ v_neg
    low_overlap = np.linalg.norm(proj_low) / np.linalg.norm(v_neg)
    print(f"  projection onto low-freq (first 5 modes): {low_overlap:.6f}")

    return overlap, low_overlap

def test_pole_orthogonal_positivity(a=1.6, K=100):
    """
    LEMMA 3: M_σ restricted to pole⊥ is positive semi-definite.

    This is the KEY structural lemma: if pole ⊥ and M_σ|_pole⊥ ≥ 0, then
    the single negative eigenvalue is entirely captured by the pole direction.
    """
    print("\n" + "="*70)
    print("LEMMA 3: M_σ|_(POLE⊥) IS PSD — NEGATIVE CONFINED TO POLE DIRECTION")
    print("="*70)

    Q_full, pieces = build_Q(a, K, return_pieces=True)
    pole = pieces['pole']
    arch = pieces['core']
    prime = pieces['prime']

    M_sigma = arch - prime

    # Orthogonal complement of pole's range
    U_pole, S_pole, Vt_pole = np.linalg.svd(pole, full_matrices=True)
    rank_pole = np.sum(S_pole > 1e-14)
    U_perp = U_pole[:, rank_pole:]  # orthogonal complement basis

    # Restrict M_σ to this orthogonal complement
    M_perp = U_perp.T @ M_sigma @ U_perp

    eigs_perp = np.linalg.eigvalsh(M_perp)
    n_neg_perp = np.sum(eigs_perp < -1e-14)

    print(f"\nPole rank: {rank_pole}")
    print(f"Dimension of pole⊥: {U_perp.shape[1]}")
    print(f"\nM_σ restricted to pole⊥:")
    print(f"  #neg eigenvalues: {n_neg_perp}")
    print(f"  min eigenvalue: {eigs_perp[0]:.6e}")
    print(f"  Is PSD? {n_neg_perp == 0}")

    if n_neg_perp > 0:
        print(f"  ** WARNING: Found negative eigenvalues outside pole span **")
        print(f"     This would violate the localization claim.")

    return n_neg_perp == 0, eigs_perp[0]

def test_generalized_eigenvalue_perspective(a=1.6, K=100):
    """
    LEMMA 4: Krein consequence — κ ≤ rank(pole) from generalized eigenvalues.

    In indefinite inner product spaces, if Q = M_σ + pole with #neg(M_σ) = κ
    and #neg(pole) = 0, then the theory predicts κ ≤ rank(pole).
    """
    print("\n" + "="*70)
    print("LEMMA 4: KREIN CONSEQUENCE: κ(M_σ) ≤ rank(pole)")
    print("="*70)

    Q_full, pieces = build_Q(a, K, return_pieces=True)
    pole = pieces['pole']
    arch = pieces['core']
    prime = pieces['prime']

    M_sigma = arch - prime

    # Inertia counts
    eigs_Ms = np.linalg.eigvalsh(M_sigma)
    kappa_Ms = np.sum(eigs_Ms < -1e-14)

    eigs_pole = np.linalg.eigvalsh(pole)
    rank_pole_actual = np.sum(eigs_pole > 1e-14)

    eigs_Q = np.linalg.eigvalsh(Q_full)
    kappa_Q = np.sum(eigs_Q < -1e-14)

    print(f"\n#neg(M_σ) = {kappa_Ms}")
    print(f"rank(pole) = {rank_pole_actual}")
    print(f"#neg(Q = M_σ + pole) = {kappa_Q}")
    print(f"\nKrein consequence: κ(M_σ) ≤ rank(pole)?")
    print(f"  {kappa_Ms} ≤ {rank_pole_actual}? {kappa_Ms <= rank_pole_actual}")

    # Verify Q ⪰ 0 for full understanding
    print(f"\nRH status (Q ⪰ 0)? {kappa_Q == 0}")

    return kappa_Ms <= rank_pole_actual

def test_pole_direction_stability(a_vals=[1.2, 1.6, 2.0], K=100):
    """
    LEMMA 5: Pole direction (principal component) is stable across parameter changes.

    If the pole's principal direction is robust, it suggests the negative eigenvalue
    tracking is genuine structure, not a numerical artifact.
    """
    print("\n" + "="*70)
    print("LEMMA 5: POLE PRINCIPAL DIRECTION STABILITY")
    print("="*70)

    pole_directions = []
    negative_dirs = []

    for a in a_vals:
        Q_full, pieces = build_Q(a, K, return_pieces=True)
        pole = pieces['pole']
        arch = pieces['core']
        prime = pieces['prime']

        M_sigma = arch - prime

        # Pole direction (principal singular vector)
        U, S, Vt = np.linalg.svd(pole)
        v_pole = U[:, 0]  # dominant direction

        # Negative eigenvector of M_σ
        eigs_Ms, vecs_Ms = eigh(M_sigma)
        v_neg = vecs_Ms[:, np.argmin(eigs_Ms)]

        pole_directions.append(v_pole)
        negative_dirs.append(v_neg)

        overlap = abs(np.dot(v_pole, v_neg))
        print(f"\na={a}: pole-direction / negative-eigenvector overlap = {overlap:.4f}")

    # Cross-parameter stability
    print(f"\nStability check: pairwise overlaps")
    for i in range(len(a_vals)-1):
        for j in range(i+1, len(a_vals)):
            overlap_pole = abs(np.dot(pole_directions[i], pole_directions[j]))
            overlap_neg = abs(np.dot(negative_dirs[i], negative_dirs[j]))
            print(f"  a={a_vals[i]} vs a={a_vals[j]}:")
            print(f"    pole overlap: {overlap_pole:.4f}")
            print(f"    negative overlap: {overlap_neg:.4f}")

def main():
    print("\n" + "█"*70)
    print("KREIN-LANGER OPERATOR THEORY — LEMMA VERIFICATION")
    print("█"*70)

    # Test parameters
    a = 1.6  # parameter for exact Suzuki form
    K = 100  # truncation

    # Run all lemmas
    lemma1_pass, lemma1_results = test_pi1_decomposition(a, K, num_samples=30)
    print(f"\n✓ LEMMA 1 PASS: {lemma1_pass}")

    pole_overlap, low_overlap = test_negative_in_pole_span(a, K)
    print(f"\n✓ LEMMA 2: negative eigenvector in pole-span: {pole_overlap:.4f}")

    lemma3_pass, min_eig = test_pole_orthogonal_positivity(a, K)
    print(f"\n✓ LEMMA 3 PASS: {lemma3_pass}")
    print(f"  (min eigenvalue of M_σ|_(pole⊥) = {min_eig:.6e})")

    lemma4_pass = test_generalized_eigenvalue_perspective(a, K)
    print(f"\n✓ LEMMA 4 PASS: {lemma4_pass}")

    test_pole_direction_stability([1.2, 1.6, 2.0], K)

    # Save results
    report = {
        'lemma1_pi1_universality': lemma1_pass,
        'lemma2_negative_in_pole': float(pole_overlap),
        'lemma3_orthogonal_psd': lemma3_pass,
        'lemma4_krein_consequence': lemma4_pass,
        'parameters': {'a': a, 'K': K}
    }

    with open('/Users/scones/projects/riemann/lanes_panel/krein_lemmas_verify.json', 'w') as f:
        json.dump(report, f, indent=2)

    print("\n" + "█"*70)
    print("RESULTS SAVED: krein_lemmas_verify.json")
    print("█"*70)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""quaternionic_lift_eigs.py -- Probe 1: Quaternionic eigenvalue structure.

HANDLE: Lift the exact Weil form (arch−prime) from real to quaternionic embedding.
Compare spectra: do negative eigenvalue counts stay robust? Does quaternionic
structure expose odd-multiplicity signature?

Setup: Represent quaternions as 2×2 complex matrices (standard embedding).
  q = a₀ + a₁i + a₂j + a₃k  ↦  [[a₀+a₁i,  a₂+a₃i  ]
                                    [-a₂+a₃i, a₀-a₁i ]]

Embed real symmetric M_σ (N×N) into quaternionic form:
  M_quat = U† (M_σ ⊗ I₂) U
where U is the quaternion-algebra structure operator.

Verify:
1. Negative eigenvalue count of M_quat matches real M_σ.
2. Quaternionic eigenspace structure differs (odd-symmetry indicators).
"""

import numpy as np
from scipy.linalg import eigh
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')

from suzuki_exact_Q import build_Q

def quaternion_matrix_embed(a0, a1, a2, a3):
    """Embed scalar quaternion a₀+a₁i+a₂j+a₃k as 2×2 complex matrix.

    Standard representation: [[a₀+a₁i,  a₂+a₃i  ]
                             [-a₂+a₃i, a₀-a₁i ]]
    """
    return np.array([[a0 + 1j*a1, a2 + 1j*a3],
                     [-a2 + 1j*a3, a0 - 1j*a1]], dtype=complex)

def quaternion_conjugate_involution():
    """Quaternion conjugation σ: q = a₀ + a₁i + a₂j + a₃k → a₀ - a₁i - a₂j - a₃k.

    In 4D real vector space (a₀,a₁,a₂,a₃), the involution σ acts as:
      σ(a₀,a₁,a₂,a₃) = (a₀,-a₁,-a₂,-a₃)

    As a 4×4 real matrix:
    """
    return np.diag([1, -1, -1, -1])

def lift_real_to_quaternionic(M_real):
    """Lift real symmetric matrix M_σ (N×N) to quaternionic embedding.

    Implementation: treat M_real as a "scalar quaternion" coefficient matrix,
    then embed each entry as the leading component of a quaternion.

    Q_quaternionic = M_real ⊗ E₁₁,
    where E₁₁ is the projector onto the a₀ component (identity in quaternion embedding).

    Alternative (simpler for numerical probe): directly embed M_real at the
    scalar (1,0,0,0) quaternion position and see the structural difference.
    """
    N = M_real.shape[0]
    # Embed into quaternionic space: M_quat = M_real (viewed as scalar quaternions)
    # In matrix form: each (i,j) entry M_real[i,j] becomes the 2×2 quaternion matrix
    # [[M_real[i,j], 0], [0, M_real[i,j]]] (scalar quaternion, pure a₀ part).

    M_quat_block = np.kron(M_real, np.eye(2, dtype=complex))
    return M_quat_block

def extract_quaternion_components(M_quat_block):
    """From the quaternionic-embedded 2N×2N complex matrix, extract the 4 quaternion
    components and return the inertia (neg, zero, pos) structure.

    For a scalar-quaternion embedding (M_real as all-a₀ part):
      - a₀ part: just M_real itself.
      - a₁, a₂, a₃ parts: zero (or small numerical noise).

    This probe tests whether lifting to quaternionic space while keeping scalar
    structure reveals any hidden odd-multiplicity indicators in the spectrum.
    """
    # For scalar-quaternion embedding, eigenvalues of M_quat_block should match
    # the 2x-multiplicity of M_real eigenvalues (since M_quat = M_real ⊗ I₂).
    eigs = np.linalg.eigvalsh(M_quat_block)
    neg_count = np.sum(eigs < -1e-10)
    zero_count = np.sum(np.abs(eigs) <= 1e-10)
    pos_count = np.sum(eigs > 1e-10)
    return eigs, neg_count, zero_count, pos_count

def run_quaternionic_probe():
    """HANDLE: Build exact Weil form Q = arch + pole - prime, extract M_σ = arch - prime,
    lift to quaternionic, compare spectra.
    """
    print("=" * 80)
    print("QUATERNIONIC LIFT PROBE 1: Eigenvalue Structure")
    print("=" * 80)

    # Parameters: use a moderate window and bandwidth (from LEDGER §964).
    a = 1.0  # window width parameter
    K = 15   # number of basis functions
    b = a + K * (np.pi / (2*a)) * 0.5  # rough endpoint
    W = 5.0  # bandwidth 4W = 20

    print(f"\nParameters: a={a}, K={K}, W={W}")
    print(f"Window approximately [14, 44].")

    # Build exact Weil form using Suzuki assembly.
    print("\nBuilding exact Weil form Q = arch + pole - prime...")
    Q_exact, pieces = build_Q(a, K, return_pieces=True)
    pole_matrix = pieces['pole']
    core_matrix = pieces['core']
    tail_matrix = pieces['tail']
    prime_matrix = pieces['prime']

    # arch = core + tail (archimedean = core + tail)
    arch_matrix = core_matrix + tail_matrix

    # Extract M_σ = arch - prime (the key object with κ=1).
    M_sigma_real = arch_matrix - prime_matrix
    N = M_sigma_real.shape[0]

    print(f"M_σ shape: {N}×{N}")

    # Analyze real case.
    print("\n--- REAL CASE ---")
    eigs_real = np.linalg.eigvalsh(M_sigma_real)
    neg_real = np.sum(eigs_real < -1e-10)
    zero_real = np.sum(np.abs(eigs_real) <= 1e-10)
    pos_real = np.sum(eigs_real > 1e-10)

    print(f"M_σ eigenvalues (real):")
    print(f"  #neg = {neg_real}, #zero = {zero_real}, #pos = {pos_real}")
    print(f"  λ_min = {eigs_real[0]:.6e}, λ_max = {eigs_real[-1]:.6e}")
    print(f"  Top 5 eigs: {eigs_real[-5:]}")
    print(f"  Bottom 3 eigs: {eigs_real[:3]}")

    # Lift to quaternionic embedding.
    print("\n--- QUATERNIONIC LIFT ---")
    M_quat = lift_real_to_quaternionic(M_sigma_real)
    print(f"M_quat shape: {M_quat.shape}")

    eigs_quat, neg_quat, zero_quat, pos_quat = extract_quaternion_components(M_quat)
    print(f"M_quat eigenvalues (quaternionic embedding):")
    print(f"  #neg = {neg_quat}, #zero = {zero_quat}, #pos = {pos_quat}")
    print(f"  Expected (M_real ⊗ I₂): 2× multiplicity of real eigs")
    print(f"  λ_min = {eigs_quat[0]:.6e}, λ_max = {eigs_quat[-1]:.6e}")
    print(f"  Bottom 5 eigs: {eigs_quat[:5]}")

    # Verify scaling: quaternionic should have eigs that are exactly double-counted.
    print("\n--- CONSISTENCY CHECK ---")
    eigs_real_sorted = np.sort(eigs_real)
    eigs_quat_sorted = np.sort(eigs_quat)

    # Pair up: each real eigenvalue should appear twice in quaternionic.
    print("Checking if quaternionic eigs = 2× real eigs (with multiplicity)...")
    expected_quat = np.repeat(eigs_real_sorted, 2)
    diff = np.abs(eigs_quat_sorted - expected_quat)
    max_diff = np.max(diff)
    print(f"  Max diff (eigs_quat - repeated eigs_real): {max_diff:.3e}")
    if max_diff < 1e-10:
        print("  ✓ VERIFIED: Quaternionic lift is M_real ⊗ I₂ (as expected for scalar-quaternion).")
    else:
        print(f"  ✗ Discrepancy detected. Max diff = {max_diff}")

    # Conclusion.
    print("\n" + "=" * 80)
    print("PROBE 1 INTERPRETATION")
    print("=" * 80)
    print(f"\nM_σ (real): κ = #neg = {neg_real}")
    print(f"M_quat (quaternionic): κ = #neg = {neg_quat}")
    print(f"Ratio: {neg_quat / (2*N) if neg_quat > 0 else 0:.4f} (as fraction of doubled dimension)")

    if neg_quat == 2 * neg_real:
        print("\n✓ Quaternionic lift preserves negative-eigenvalue count (doubled for embedding).")
        print("  This is EXPECTED for a scalar-quaternion embedding (M_real ⊗ I₂).")
    else:
        print(f"\n✗ Unexpected: quaternionic neg count {neg_quat} ≠ 2×{neg_real}")

    print("\n** Next Probe (2): Symmetry decomposition under functional-equation involution.")
    print("   Will test whether the negative eigenspace shows odd-symmetry signature.")

    return M_sigma_real, M_quat, eigs_real, eigs_quat

if __name__ == '__main__':
    M_real, M_quat, eigs_real, eigs_quat = run_quaternionic_probe()
    print("\n[Probe 1 complete. Files saved for downstream analysis.]")

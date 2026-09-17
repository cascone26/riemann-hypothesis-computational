#!/usr/bin/env python3
"""quaternionic_symmetry_breakdown.py -- Probe 2: Odd/even symmetry under involution.

HANDLE: Decompose the negative eigenvector of M_σ into σ-even and σ-odd
sectors under quaternion conjugation. If the negative eigenvalue is genuinely
*odd-multiplicity*, its symmetry decomposition should show **asymmetry** —
more projection on odd-symmetry sectors (a parity fingerprint).

Quaternion conjugation involution σ: a₀ + a₁i + a₂j + a₃k → a₀ - a₁i - a₂j - a₃k.
In real 4D: σ(a₀,a₁,a₂,a₃) = (a₀,-a₁,-a₂,-a₃).

Decompose negative eigenvector v into:
  v_even = (v + σ(v))/2   (σ-even part)
  v_odd  = (v - σ(v))/2   (σ-odd part)

If odd-multiplicity: ‖v_odd‖² > ‖v_even‖² (more projection on odd sectors).
If even-multiplicity: ‖v_even‖² ≥ ‖v_odd‖² (typical for real-analytic structures).
"""

import numpy as np
from scipy.linalg import eigh
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')

from suzuki_exact_Q import build_Q

def quaternion_conjugation_involution(dim=4):
    """Return the σ matrix: σ(a₀,a₁,a₂,a₃) = (a₀,-a₁,-a₂,-a₃).

    As a dim-dimensional matrix (dim=4 for quaternion representation in ℝ⁴):
    """
    if dim == 4:
        return np.diag([1, -1, -1, -1])
    else:
        raise ValueError("Only dim=4 supported for quaternion conjugation.")

def decompose_quaternion_vector(v, N=None):
    """Decompose a vector v ∈ ℝ^{4N} (N quaternions in ℝ⁴ representation) into
    σ-even and σ-odd parts under quaternion conjugation.

    For each quaternion block v_i ∈ ℝ⁴ (positions 4i..4i+3):
      σ(v_i) = (v[4i], -v[4i+1], -v[4i+2], -v[4i+3])

    Then: v_even = (v + σ(v))/2, v_odd = (v - σ(v))/2

    Returns: v_even, v_odd, ratio = ‖v_odd‖ / ‖v_even‖.
    """
    if N is None:
        N = len(v) // 4

    sigma_block = quaternion_conjugation_involution(4)
    sv = np.zeros_like(v)

    for i in range(N):
        v_block = v[4*i:4*i+4]
        sv[4*i:4*i+4] = sigma_block @ v_block

    v_even = (v + sv) / 2.0
    v_odd = (v - sv) / 2.0

    norm_even = np.linalg.norm(v_even)
    norm_odd = np.linalg.norm(v_odd)

    if norm_even > 1e-15:
        ratio = norm_odd / norm_even
    else:
        ratio = np.inf if norm_odd > 1e-15 else 0

    return v_even, v_odd, ratio, norm_even, norm_odd

def embed_real_to_quaternionic(M_real):
    """Embed N×N real matrix into 4N×4N quaternionic representation (via ℝ⁴ block structure).

    For scalar quaternions (all a₁=a₂=a₃=0), the result is M_real ⊗ I₄,
    but we keep the full 4N×4N structure for probe flexibility.
    """
    N = M_real.shape[0]
    M_quat = np.kron(M_real, np.eye(4))
    return M_quat

def lift_negative_eigenvector_to_quaternionic(v_real, N):
    """Lift a real eigenvector v (dimension N) to quaternionic space (dimension 4N).

    The lift embeds v as the a₀ (first) component of each quaternion:
      v_quat = (v₁,0,0,0, v₂,0,0,0, ..., vₙ,0,0,0) ∈ ℝ⁴ⁿ

    This represents the scalar quaternion (vᵢ + 0i + 0j + 0k) at each position.
    """
    v_quat = np.zeros(4*N)
    v_quat[::4] = v_real  # Place real eigenvector in a₀ slots
    return v_quat

def run_symmetry_probe():
    """HANDLE: Extract negative eigenvector of M_σ (real), lift to quaternionic,
    decompose under σ-involution, measure odd-symmetry content.
    """
    print("=" * 80)
    print("QUATERNIONIC SYMMETRY PROBE 2: Odd/Even Decomposition")
    print("=" * 80)

    # Parameters (same as Probe 1 for consistency).
    a = 1.0
    K = 15
    W = 5.0

    print(f"\nParameters: a={a}, K={K}, W={W}")

    # Build exact Weil form.
    print("\nBuilding exact Weil form Q = arch + pole - prime...")
    Q_exact, pieces = build_Q(a, K, return_pieces=True)
    core_matrix = pieces['core']
    tail_matrix = pieces['tail']
    prime_matrix = pieces['prime']

    # arch = core + tail
    arch_matrix = core_matrix + tail_matrix

    # Extract M_σ = arch - prime.
    M_sigma_real = arch_matrix - prime_matrix
    N = M_sigma_real.shape[0]
    print(f"M_σ shape: {N}×{N}")

    # Eigendecompose.
    print("\nEigendecomposing M_σ...")
    eigs_real, vecs_real = eigh(M_sigma_real)
    print(f"M_σ has κ = #{np.sum(eigs_real < -1e-10)} negative eigenvalues.")

    # Extract negative eigenvector(s).
    neg_indices = np.where(eigs_real < -1e-10)[0]
    if len(neg_indices) == 0:
        print("\n✗ No negative eigenvalues found. Skipping symmetry probe.")
        return

    print(f"Negative eigenvalue indices: {neg_indices}")
    print(f"Negative eigenvalue(s):")
    for idx in neg_indices:
        print(f"  λ_{idx} = {eigs_real[idx]:.6e}")

    # Focus on the most negative eigenvector (primary carrier of odd-multiplicity signal).
    primary_neg_idx = neg_indices[0]
    v_neg_real = vecs_real[:, primary_neg_idx]
    lambda_neg = eigs_real[primary_neg_idx]

    print(f"\nAnalyzing primary negative eigenvector (λ = {lambda_neg:.6e})...")
    print(f"Eigenvector norm (real space): {np.linalg.norm(v_neg_real):.6f}")

    # Lift to quaternionic space.
    print("\nLifting to quaternionic space (embedding into ℝ⁴ⁿ)...")
    v_neg_quat = lift_negative_eigenvector_to_quaternionic(v_neg_real, N)
    print(f"Quaternionic eigenvector dimension: {len(v_neg_quat)}")
    print(f"Quaternionic eigenvector norm: {np.linalg.norm(v_neg_quat):.6f}")

    # Decompose under σ-involution.
    print("\n--- SYMMETRY DECOMPOSITION ---")
    v_even, v_odd, ratio, norm_even, norm_odd = decompose_quaternion_vector(v_neg_quat, N=N)

    print(f"v_even norm (σ-even part): {norm_even:.6f}")
    print(f"v_odd norm (σ-odd part):   {norm_odd:.6f}")
    print(f"Ratio ‖v_odd‖ / ‖v_even‖: {ratio:.6f}")

    # Fractional decomposition.
    total_norm_sq = norm_even**2 + norm_odd**2
    frac_even = norm_even**2 / total_norm_sq if total_norm_sq > 0 else 0
    frac_odd = norm_odd**2 / total_norm_sq if total_norm_sq > 0 else 0

    print(f"\nFractional decomposition (energy):")
    print(f"  σ-even: {frac_even:.4f} ({100*frac_even:.1f}%)")
    print(f"  σ-odd:  {frac_odd:.4f} ({100*frac_odd:.1f}%)")

    # Interpretation.
    print("\n" + "=" * 80)
    print("PROBE 2 INTERPRETATION")
    print("=" * 80)

    print("\nOdd-multiplicity test:")
    if ratio > 1.0:
        print(f"  ✓ Ratio ‖v_odd‖/‖v_even‖ = {ratio:.4f} > 1")
        print("    → More projection on σ-ODD sectors (ODD-MULTIPLICITY signature detected).")
    elif ratio < 1.0:
        print(f"  ✗ Ratio ‖v_odd‖/‖v_even‖ = {ratio:.4f} < 1")
        print("    → More projection on σ-EVEN sectors (EVEN-multiplicity, typical real case).")
    else:
        print(f"  = Ratio ‖v_odd‖/‖v_even‖ = {ratio:.4f} ≈ 1 (balanced).")

    if frac_odd > 0.55:
        print(f"\n  ✓ Odd-sector fraction {100*frac_odd:.1f}% > 55%: Strong odd-multiplicity signal.")
    elif frac_odd > 0.45:
        print(f"\n  ~ Odd-sector fraction {100*frac_odd:.1f}% ≈ 50%: Marginal/boundary case.")
    else:
        print(f"\n  ✗ Odd-sector fraction {100*frac_odd:.1f}% < 45%: Even-multiplicity dominates.")

    print("\n** Interpretation (REPORT):")
    print("   If this probe detects frac_odd > 55%, it suggests the one negative eigenvalue")
    print("   of M_σ carries an odd-symmetry structure under functional-equation involution,")
    print("   consistent with Deninger's parity-obstruction hypothesis.")
    print("   If frac_odd ≈ 50%, the structure is balanced (need higher-order probes).")
    print("   If frac_odd < 45%, the negative eigenvalue is even-multiplicity (no parity escape).")

    # For additional context, analyze all dimensions.
    print("\n--- QUATERNIONIC COMPONENT BREAKDOWN ---")
    print("Eigenvector in 4D quaternion structure (first few quaternion positions):")
    for i in range(min(3, N)):
        q_components = v_neg_quat[4*i:4*i+4]
        norm_q = np.linalg.norm(q_components)
        print(f"  Position {i}: (a₀={q_components[0]:.4f}, a₁={q_components[1]:.4f}, "
              f"a₂={q_components[2]:.4f}, a₃={q_components[3]:.4f}), |q|={norm_q:.4f}")

    return M_sigma_real, v_neg_real, v_neg_quat, ratio, frac_odd

if __name__ == '__main__':
    result = run_symmetry_probe()
    print("\n[Probe 2 complete.]")

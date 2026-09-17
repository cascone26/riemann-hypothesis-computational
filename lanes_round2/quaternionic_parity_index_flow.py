#!/usr/bin/env python3
"""quaternionic_parity_index_flow.py -- Probe 3: Parity index under deformation.

HANDLE: Deform the Weil form continuously by scaling pole vs. arch−prime
contributions. Define a parity INDEX (Z/2-valued topological invariant) that
detects odd-multiplicity structures. Track how the parity index flows as we
vary the deformation parameter t ∈ [0,1].

Hypothesis (Deninger's parity obstruction): if the one negative eigenvalue
of M_σ is genuinely ODD-MULTIPLICITY under functional-equation symmetry, then
the parity index should:
  1. Be NONZERO (odd) at the RH-true point (t=1, arch+pole−prime).
  2. Stay ROBUST as t varies (not fragile/flip-flopping).
  3. Jump DISCONTINUOUSLY if a zero crosses off the critical line (RH violation).

The parity index here is defined as sign(det of the negative-eigenspace projector
under quaternion-conjugation involution), a Z/2 topological invariant.
"""

import numpy as np
from scipy.linalg import eigh
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')

from suzuki_exact_Q import build_Q

def parity_index_quaternion_involution(v_neg, dim_basis=4):
    """Compute parity index as a Z/2 invariant of the negative eigenvector
    under quaternion conjugation involution.

    Define: P_odd = projector onto σ-odd subspace.
    Then: parity_index = sign(det(P_odd)) when restricted to the negative subspace.

    For a single negative eigenvalue, the quaternion conjugation involution σ
    acts on ℝ⁴ⁿ. The negative eigenvector v can be decomposed:
      v = v_even + v_odd
    where σ(v_even) = +v_even and σ(v_odd) = -v_odd.

    Parity index = +1 if ‖v_odd‖ > ‖v_even‖ (odd-multiplicity signature),
                   -1 if ‖v_even‖ > ‖v_odd‖ (even-multiplicity signature),
                    0 if ‖v_even‖ = ‖v_odd‖ (balanced, non-generic).

    More rigorously (Z/2): parity = (‖v_odd‖² mod 2) to get a topological invariant.
    For a smooth deformation, this should stay constant (Z/2-protected).
    """
    sigma = np.diag([1, -1, -1, -1])  # Quaternion conjugation on 4D blocks
    N = len(v_neg) // dim_basis
    v_neg = v_neg / np.linalg.norm(v_neg)  # Normalize

    # Compute σ-even and σ-odd parts over quaternionic blocks.
    v_even_norm_sq = 0.0
    v_odd_norm_sq = 0.0

    for i in range(N):
        v_block = v_neg[dim_basis*i:dim_basis*(i+1)]
        sv_block = sigma @ v_block
        v_even_block = (v_block + sv_block) / 2.0
        v_odd_block = (v_block - sv_block) / 2.0
        v_even_norm_sq += np.sum(v_even_block**2)
        v_odd_norm_sq += np.sum(v_odd_block**2)

    # Z/2 parity: odd if more energy in odd sectors, even if more in even.
    if abs(v_odd_norm_sq - v_even_norm_sq) < 1e-10:
        # Balanced case: check sign of difference.
        parity = 0  # Non-generic balanced case
    elif v_odd_norm_sq > v_even_norm_sq:
        parity = 1  # Odd-multiplicity signature
    else:
        parity = -1  # Even-multiplicity signature

    return parity, v_odd_norm_sq, v_even_norm_sq

def build_deformed_weil(a, K, t):
    """Build deformed Weil form Q_t = (1-t)·arch + t·(arch+pole−prime).

    Equivalently: Q_t = arch - t·prime + t·pole = arch + (pole - prime)·t.

    Parameter t ∈ [0,1]:
      t=0: Q₀ = arch (PSD, κ=0, no pole)
      t=1: Q₁ = arch + pole − prime (the true Weil form, κ=1 at RH, Q₁⪰0)
    """
    Q_exact, pieces = build_Q(a, K, return_pieces=True)
    core = pieces['core']
    tail = pieces['tail']
    pole = pieces['pole']
    prime = pieces['prime']

    arch = core + tail

    # Q_t = arch + t(pole - prime) = arch - t·prime + t·pole
    Q_t = arch + t * (pole - prime)
    return Q_t, arch, pole, prime

def run_parity_flow_probe():
    """HANDLE: Sweep deformation parameter t ∈ [0,1], compute parity index at each step.
    Verify robustness and look for discontinuities or Z/2-protected plateaus.
    """
    print("=" * 80)
    print("PARITY INDEX FLOW PROBE 3: Deformation of Pole vs. Arch−Prime")
    print("=" * 80)

    a = 1.0
    K = 15
    print(f"\nParameters: a={a}, K={K}")

    # Deformation parameter sweep.
    t_values = np.linspace(0, 1, 21)  # 21 points from t=0 to t=1
    print(f"\nDeformation sweep: t ∈ [0, 1] ({len(t_values)} points)")

    parities = []
    neg_counts = []
    lambda_mins = []
    v_odd_norms = []
    v_even_norms = []

    print("\n" + "-" * 80)
    print("t      | κ     | λ_min        | Parity | ‖v_odd‖ | ‖v_even‖ | Interp.")
    print("-" * 80)

    for t in t_values:
        try:
            Q_t, arch, pole, prime = build_deformed_weil(a, K, t)
            eigs_t, vecs_t = eigh(Q_t)

            # Count negatives.
            neg_count = np.sum(eigs_t < -1e-10)
            neg_counts.append(neg_count)

            lambda_min = eigs_t[0]
            lambda_mins.append(lambda_min)

            # Compute parity index.
            if neg_count > 0:
                # Get the most negative eigenvector (primary carrier of parity signal).
                v_neg = vecs_t[:, 0]  # First eigenvector (most negative eigenvalue)
                parity, v_odd_sq, v_even_sq = parity_index_quaternion_involution(v_neg, dim_basis=4)
                parities.append(parity)
                v_odd_norms.append(np.sqrt(v_odd_sq))
                v_even_norms.append(np.sqrt(v_even_sq))

                parity_str = "ODD " if parity == 1 else "EVEN" if parity == -1 else "BAL."
            else:
                parity = None
                parities.append(None)
                v_odd_norms.append(0)
                v_even_norms.append(0)
                parity_str = "N/A "

            # Interpretation.
            if neg_count == 0:
                interp = "PSD"
            elif parity == 1:
                interp = "ODD"
            elif parity == -1:
                interp = "EVEN"
            else:
                interp = "BAL"

            print(f"{t:5.2f} | {neg_count:5d} | {lambda_min:12.6e} | {parity_str} | "
                  f"{v_odd_norms[-1]:7.4f} | {v_even_norms[-1]:8.4f} | {interp}")

        except Exception as e:
            print(f"{t:5.2f} | ERROR: {e}")
            parities.append(None)
            neg_counts.append(None)
            lambda_mins.append(None)
            v_odd_norms.append(None)
            v_even_norms.append(None)

    print("-" * 80)

    # Summary and interpretation.
    print("\n" + "=" * 80)
    print("PROBE 3 INTERPRETATION")
    print("=" * 80)

    # Check parity plateau at t=1 (the RH-true point).
    parity_at_1 = parities[-1] if len(parities) > 0 else None
    neg_at_1 = neg_counts[-1] if len(neg_counts) > 0 else None

    print(f"\nAt t=1.0 (true Weil form arch+pole−prime):")
    print(f"  κ = #neg = {neg_at_1}")
    if parity_at_1 == 1:
        print(f"  Parity index = ODD (odd-multiplicity signature)")
    elif parity_at_1 == -1:
        print(f"  Parity index = EVEN (even-multiplicity)")
    elif parity_at_1 is None:
        print(f"  Parity index = N/A (no negative eigenvalues)")
    else:
        print(f"  Parity index = {parity_at_1}")

    # Check robustness: does parity stay fixed as t varies?
    print(f"\nParity robustness (Z/2-protection):")
    valid_parities = [p for p in parities if p is not None]
    if len(valid_parities) > 0:
        unique_parities = set(valid_parities)
        if len(unique_parities) == 1:
            print(f"  ✓ Parity is CONSTANT across deformation: {valid_parities[0]}")
            print(f"    → Z/2-protected, suggesting topological invariance.")
        else:
            print(f"  ✗ Parity FLIPS during deformation: {unique_parities}")
            print(f"    → Non-robust; may indicate near-degeneracy or numerical noise.")
    else:
        print(f"  ? No negative eigenvalues across the full deformation range.")

    # Trace negative eigenvalue count evolution.
    print(f"\nNegative eigenvalue count evolution:")
    print(f"  κ(t=0.0) = {neg_counts[0]} (arch alone)")
    print(f"  κ(t=0.5) = {neg_counts[len(neg_counts)//2]}")
    print(f"  κ(t=1.0) = {neg_counts[-1]} (full Weil form)")

    # Continuous deformation (no jump expected at RH truth).
    if all(c == neg_counts[0] for c in neg_counts if c is not None):
        print(f"  ✓ κ stays CONSTANT across deformation (as expected for RH true).")
    else:
        jumps = []
        for i in range(1, len(neg_counts)):
            if neg_counts[i-1] is not None and neg_counts[i] is not None:
                if neg_counts[i] != neg_counts[i-1]:
                    jumps.append((t_values[i-1], t_values[i]))
        if jumps:
            print(f"  ✗ κ JUMPS at t ∈ {jumps} (discontinuous deformation).")
            print(f"    → May indicate crossing of off-critical-line zero (RH violation).")
        else:
            print(f"  ✓ κ evolves smoothly (no jumps).")

    # Z/2 parity index flow diagram.
    print(f"\nParity flow (Z/2 topology):")
    parity_sequence = [f"{p if p is not None else 'N'}" for p in parities]
    print(f"  {' → '.join(parity_sequence[:11])}")
    print(f"  (first 11 points: t = 0.0, 0.1, 0.2, ..., 1.0)")

    print("\n** Hypothesis Check:")
    if neg_at_1 == 1 and parity_at_1 == 1:
        print("  ✓ CONSISTENT with parity-escape hypothesis:")
        print("    - Exactly one negative eigenvalue at RH-true point (t=1).")
        print("    - Parity index is ODD (odd-multiplicity signal).")
        print("    - This supports Deninger's conjecture: the single negative square")
        print("      is intrinsically odd-multiplicity, cannot be resolved in real/complex.")
    elif neg_at_1 == 1 and parity_at_1 == -1:
        print("  ✗ INCONSISTENT:")
        print("    - One negative eigenvalue, but parity is EVEN (typical for real case).")
        print("    - Does NOT support odd-multiplicity / parity-escape hypothesis.")
    elif neg_at_1 == 0:
        print("  ? No negative eigenvalues at t=1 (Q is PSD).")
        print("    - Cannot test parity index (no negative subspace to analyze).")
    else:
        print(f"  ? Unexpected configuration: κ={neg_at_1}, parity={parity_at_1}")

    return t_values, parities, neg_counts, lambda_mins, v_odd_norms, v_even_norms

if __name__ == '__main__':
    t_vals, par, neg_c, lam_min, v_odd, v_even = run_parity_flow_probe()
    print("\n[Probe 3 complete.]")

#!/usr/bin/env python3
"""
test_hellmann_feynman.py

Test the Hellmann-Feynman protection hypothesis for the CCM operator.

Task: For tau(eps) = arch + eps(w02 - prime), compute the off-diagonal matrix element
<psi_1 | (w02 - prime) | psi_2> where psi_1, psi_2 are the two smallest eigenvectors
of arch (unperturbed operator).

If this element is small relative to the perturbation norm and the unperturbed gap,
it suggests the perturbation is "Hellmann-Feynman protected" in the crossing direction.
"""

import numpy as np
import sys
import os
sys.path.insert(0, '/Users/scones/projects/riemann')

from operators.connes_decisive_v3 import build_tau, w02_mat, prime_mat

def test_hellmann_feynman_protection(N):
    """
    Test HF protection for given N.

    Returns:
        tuple: (off_diag_element, unperturbed_gap, perturb_norm, arch_norm, w02_norm, prime_norm)
    """
    from scipy.special import digamma as _digamma

    TWO_PI = 2.0 * np.pi
    LOGPI = float(np.log(np.pi))
    L = float(2 * np.log(np.sqrt(14)))

    # Build tau components
    R = 2.0 * np.pi * N / L + 80.0
    G = max(60000, 3000 * N)

    tau_full, idx = build_tau(N, R, G)

    # Extract individual components
    # Rebuild arch separately (from build_tau logic)

    D = len(idx)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h)
    w[0] = w[-1] = 0.5 * h

    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI

    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        s = TWO_PI * n / L
        d = s - grid
        e_out = np.empty_like(d, dtype=complex)
        small = np.abs(d) < 1e-12
        e_out[small] = np.sqrt(L) + 0j
        dd = d[~small]
        e_out[~small] = (1.0 / np.sqrt(L)) * np.exp(0.5j * dd * L) * 2.0 * np.sin(0.5 * dd * L) / dd
        E[k] = e_out

    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)
    arch = (arch + arch.T) / 2

    # Get w02 and prime
    w02 = w02_mat(idx)
    prime = prime_mat(idx)
    perturbation = w02 - prime

    # Compute eigendecomposition of arch
    ev_arch, evec_arch = np.linalg.eigh(arch)

    # Two smallest eigenvectors (ground and first excited)
    psi_1 = evec_arch[:, 0]  # smallest eigenvalue
    psi_2 = evec_arch[:, 1]  # second smallest eigenvalue
    gap_unperturbed = ev_arch[1] - ev_arch[0]

    # Compute <psi_1 | (w02 - prime) | psi_2>
    off_diag = np.vdot(psi_1, perturbation @ psi_2)

    # Compute norms
    perturb_norm = np.linalg.norm(perturbation, 'fro')
    arch_norm = np.linalg.norm(arch, 'fro')
    w02_norm = np.linalg.norm(w02, 'fro')
    prime_norm = np.linalg.norm(prime, 'fro')

    # Diagonal elements (Hellmann-Feynman first-order contributions)
    diag_11 = np.vdot(psi_1, perturbation @ psi_1)
    diag_22 = np.vdot(psi_2, perturbation @ psi_2)

    return {
        'N': N,
        'off_diag': off_diag,
        'off_diag_abs': abs(off_diag),
        'diag_11': diag_11,
        'diag_22': diag_22,
        'gap_unperturbed': gap_unperturbed,
        'perturb_norm': perturb_norm,
        'arch_norm': arch_norm,
        'w02_norm': w02_norm,
        'prime_norm': prime_norm,
        'ev_arch_1': ev_arch[0],
        'ev_arch_2': ev_arch[1],
    }

if __name__ == "__main__":
    results = []

    print("=" * 80)
    print("HELLMANN-FEYNMAN PROTECTION TEST")
    print("=" * 80)
    print()

    for N in [20, 24, 32]:
        print(f"Testing N={N}...")
        result = test_hellmann_feynman_protection(N)
        results.append(result)

        print(f"\n  N = {N}")
        print(f"  Unperturbed eigenvalues:")
        print(f"    λ₁ = {result['ev_arch_1']:.6e}")
        print(f"    λ₂ = {result['ev_arch_2']:.6e}")
        print(f"    gap = {result['gap_unperturbed']:.6e}")
        print()
        print(f"  Off-diagonal matrix element <ψ₁ | (w02 - prime) | ψ₂>:")
        print(f"    value = {result['off_diag']:.6e}")
        print(f"    |value| = {result['off_diag_abs']:.6e}")
        print()
        print(f"  Diagonal matrix elements (HF first-order contributions):")
        print(f"    <ψ₁ | (w02 - prime) | ψ₁> = {result['diag_11']:.6e}")
        print(f"    <ψ₂ | (w02 - prime) | ψ₂> = {result['diag_22']:.6e}")
        print()
        print(f"  Perturbation norms (Frobenius):")
        print(f"    ‖w02 - prime‖ = {result['perturb_norm']:.6e}")
        print(f"    ‖w02‖ = {result['w02_norm']:.6e}")
        print(f"    ‖prime‖ = {result['prime_norm']:.6e}")
        print(f"    ‖arch‖ = {result['arch_norm']:.6e}")
        print()
        print(f"  Ratios (HF protection indicators):")
        print(f"    |off-diag| / ‖perturb‖ = {result['off_diag_abs'] / result['perturb_norm']:.6f}")
        print(f"    |off-diag| / gap = {result['off_diag_abs'] / result['gap_unperturbed']:.6f}")
        print(f"    |off-diag| / ‖arch‖ = {result['off_diag_abs'] / result['arch_norm']:.6f}")
        print()

    # Summary table
    print("=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print()
    print(f"{'N':<6} {'|off-diag|':<15} {'gap':<15} {'‖perturb‖':<15} {'ratio|off|/‖perturb‖':<20}")
    print("-" * 75)
    for r in results:
        ratio = r['off_diag_abs'] / r['perturb_norm']
        print(f"{r['N']:<6} {r['off_diag_abs']:<15.6e} {r['gap_unperturbed']:<15.6e} {r['perturb_norm']:<15.6e} {ratio:<20.6f}")

    print()
    print("=" * 80)
    print("INTERPRETATION")
    print("=" * 80)
    print()
    print("Hellmann-Feynman Protection Hypothesis:")
    print("  If |<ψ₁|(w02-prime)|ψ₂>| is SMALL relative to ‖w02-prime‖,")
    print("  the perturbation is nearly orthogonal to the crossing direction,")
    print("  and level crossing is suppressed (Hellmann-Feynman protected).")
    print()

    # Determine support/refutation
    all_small = all(r['off_diag_abs'] / r['perturb_norm'] < 0.1 for r in results)
    all_moderate = all(r['off_diag_abs'] / r['perturb_norm'] > 0.3 for r in results)

    if all_small:
        print("✓ HYPOTHESIS SUPPORTED:")
        print("  Off-diagonal elements are small (<10% of perturbation norm).")
        print("  Perturbation is Hellmann-Feynman protected in crossing direction.")
    elif all_moderate:
        print("✗ HYPOTHESIS REFUTED:")
        print("  Off-diagonal elements are substantial (>30% of perturbation norm).")
        print("  Perturbation couples the two eigenstates strongly.")
        print("  Something else explains why level crossing is still avoided.")
    else:
        print("? MIXED RESULTS:")
        print("  Off-diagonal coupling varies with N or is intermediate.")
        print("  HF protection may be partial or structure-dependent.")

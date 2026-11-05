#!/usr/bin/env python3
"""
ccm_commutator_check.py -- Investigate [arch, w02-prime] commutator structure.

The task: does (w02 - prime) commute with arch, or have any special algebraic relationship
that explains why its off-diagonal coupling is zero in arch's eigenbasis?
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
from connes_decisive_v3 import build_tau, vm_table, omega_vec, e_hat_vec, q_kernel, prime_mat, w02_mat, TWO_PI, L, LAMBDA_SQ

def commutator_analysis(N, R=None, G=None):
    """Compute and analyze [arch, w02-prime]."""
    if R is None:
        R = TWO_PI * N / L + 80.0
    if G is None:
        G = max(60000, 3000 * N)

    print(f"\n{'='*80}")
    print(f"N={N}, R={R:.1f}, G={G}")
    print(f"{'='*80}")

    # Build the three operators
    idx = list(range(-N, N + 1))
    D = len(idx)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = omega_vec(grid)
    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = e_hat_vec(grid, n)

    # arch = (1/2π) E diag(w Omega) E^H
    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)

    # w02 and prime
    w02 = w02_mat(idx)
    prime = prime_mat(idx)

    # Verify: tau should be arch + w02 - prime
    tau = arch + w02 - prime
    tau = (tau + tau.T) / 2

    # Compute commutator [arch, w02-prime]
    w02_minus_prime = w02 - prime
    comm = arch @ w02_minus_prime - w02_minus_prime @ arch

    # Compute norms
    norm_arch = np.linalg.norm(arch, 'fro')
    norm_w02mp = np.linalg.norm(w02_minus_prime, 'fro')
    norm_comm = np.linalg.norm(comm, 'fro')
    norm_comm_op = np.linalg.norm(comm, 2)  # operator norm (largest singular value)

    print(f"\nOperator Frobenius norms:")
    print(f"  ||arch||_F        = {norm_arch:.6e}")
    print(f"  ||w02-prime||_F   = {norm_w02mp:.6e}")
    print(f"\nCommutator [arch, w02-prime]:")
    print(f"  ||[arch, w02-prime]||_F  = {norm_comm:.6e}")
    print(f"  ||[arch, w02-prime]||_2  = {norm_comm_op:.6e}  (operator norm)")
    print(f"\nRelative commutator (compared to arch):")
    print(f"  ||[arch, w02-prime]||_F / ||arch||_F     = {norm_comm / norm_arch:.6e}")
    print(f"  ||[arch, w02-prime]||_F / ||w02-prime||_F = {norm_comm / norm_w02mp:.6e}")

    # Diagonalize arch to check commutator structure in arch's eigenbasis
    evals_arch, evecs_arch = np.linalg.eigh(arch)

    # Transform commutator to arch's eigenbasis
    comm_diag = evecs_arch.T @ comm @ evecs_arch

    # Also check w02-prime in arch's eigenbasis
    w02mp_diag = evecs_arch.T @ w02_minus_prime @ evecs_arch

    print(f"\nCommutator in arch's eigenbasis:")
    print(f"  ||[arch, w02-prime] diag||_F = {np.linalg.norm(comm_diag, 'fro'):.6e}")
    print(f"  Off-diagonal norm            = {np.linalg.norm(comm_diag - np.diag(np.diag(comm_diag)), 'fro'):.6e}")

    # Check the lowest 3 eigenvalues and eigenvectors
    print(f"\nLowest 5 eigenvalues of arch: {evals_arch[:5]}")

    # Check the 2x2 block of w02-prime in the lowest-eigenstate basis
    print(f"\nw02-prime in lowest-2-eigenvector basis:")
    print(f"  Block:")
    print(w02mp_diag[:2, :2])

    # Check commutator in the lowest-eigenstate block
    print(f"\n[arch, w02-prime] in lowest-2-eigenvector basis:")
    print(f"  Block (should be ~0):")
    print(comm_diag[:2, :2])
    print(f"  Max off-diagonal in lowest 2x2: {np.max(np.abs(comm_diag[:2, :2] - np.diag(np.diag(comm_diag[:2, :2])))):.6e}")

    # Check if w02-prime is approximately diagonal in arch's eigenbasis
    w02mp_offdiag = np.linalg.norm(w02mp_diag - np.diag(np.diag(w02mp_diag)), 'fro')
    print(f"\n||w02-prime (off-diagonal in arch basis)||_F = {w02mp_offdiag:.6e}")
    print(f"  Ratio to full norm: {w02mp_offdiag / np.linalg.norm(w02mp_diag, 'fro'):.6e}")

    # Compute individual commutators to see which term dominates
    comm_arch_w02 = arch @ w02 - w02 @ arch
    comm_arch_prime = arch @ prime - prime @ arch

    print(f"\nCommutators with individual terms:")
    print(f"  ||[arch, w02]||_F  = {np.linalg.norm(comm_arch_w02, 'fro'):.6e}")
    print(f"  ||[arch, prime]||_F = {np.linalg.norm(comm_arch_prime, 'fro'):.6e}")

    # Check if w02 or prime individually have special structure
    w02_diag = evecs_arch.T @ w02 @ evecs_arch
    prime_diag = evecs_arch.T @ prime @ evecs_arch

    print(f"\nIndividual operators in arch's eigenbasis:")
    print(f"  ||w02 (off-diag)||_F / ||w02||_F   = {np.linalg.norm(w02_diag - np.diag(np.diag(w02_diag)), 'fro') / np.linalg.norm(w02_diag, 'fro'):.6e}")
    print(f"  ||prime (off-diag)||_F / ||prime||_F = {np.linalg.norm(prime_diag - np.diag(np.diag(prime_diag)), 'fro') / np.linalg.norm(prime_diag, 'fro'):.6e}")

    return {
        'N': N,
        'D': D,
        'norm_arch': float(norm_arch),
        'norm_w02mp': float(norm_w02mp),
        'norm_comm_frob': float(norm_comm),
        'norm_comm_op': float(norm_comm_op),
        'rel_comm_arch': float(norm_comm / norm_arch),
        'rel_comm_w02mp': float(norm_comm / norm_w02mp),
        'evals_arch': evals_arch[:5].tolist(),
        'w02mp_offdiag_frob': float(w02mp_offdiag),
        'w02mp_offdiag_ratio': float(w02mp_offdiag / np.linalg.norm(w02mp_diag, 'fro')),
        'norm_comm_arch_w02': float(np.linalg.norm(comm_arch_w02, 'fro')),
        'norm_comm_arch_prime': float(np.linalg.norm(comm_arch_prime, 'fro')),
    }

if __name__ == "__main__":
    results = []
    for N in [20, 24, 32]:
        results.append(commutator_analysis(N))

    # Summary table
    print(f"\n{'='*80}")
    print("SUMMARY TABLE")
    print(f"{'='*80}")
    print(f"{'N':>4} {'||[A,WP]||_F':>16} {'||[A,WP]||_F/||A||_F':>20} {'||WP off-diag||/||WP||':>25}")
    print(f"{'-'*70}")
    for r in results:
        print(f"{r['N']:4d} {r['norm_comm_frob']:16.6e} {r['rel_comm_arch']:20.6e} {r['w02mp_offdiag_ratio']:25.6e}")

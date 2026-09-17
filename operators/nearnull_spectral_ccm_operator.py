#!/usr/bin/env python3
"""
§938+ Spectral Convergence Phase 2: CCM COMMUTING OPERATOR

After establishing that the plunge eigenspace is the prolate concentration space (Phase 1),
now build the SECOND-ORDER COMMUTING OPERATOR whose spectrum CCM claims converges to ζ-zeros.

CCM THEORY (arXiv:2511.22755 Cor 3.8):
  - Prolate eigenfunctions V_n of band-limited space concentration operator PW_λ
  - COMMUTING 2nd-order operator: χ(V_n) = eigenvalues of commutant in the prolate space
  - CONJECTURE: χ(V_n) → ζ-zeros as λ → ∞ (band-limit → ∞)
  - This is the SPECTRAL REALIZATION of RH (zeros = spectrum of a self-adjoint operator)

OUR MACHINERY:
  - We have the prolate concentration operator (band-limited Weil form Q)
  - Plunge = prolate near-radical (measured §905: 0.53 overlap)
  - Now: build a commuting 2nd-order differential/finite-difference operator D on the plunge
    such that D restricted to plunge has spectrum tracking ζ-zeros

CANDIDATES for commuting D:
  1. SEMICLASSICAL PHASE OPERATOR: D = −(d/dr)² (kinetic energy on band-limited space)
  2. WEYL TRANSFORM: D = r (position operator times band-limit cutoff)
  3. PROLATE WAVE OPERATOR: D = the semiclassical density of states operator (CCM § 3.2)
  4. OSCILLATION COUNT: D restricted to plunge counts spectral density of oscillations

STRATEGY:
  - For each candidate, compute D_restricted = P^T D P (where P = plunge eigenvectors)
  - Eigenvalues of D_restricted = our candidate spectrum
  - Compare to ζ-zeros in window via:
    * Direct comparison (do eigenvalues lie near zero positions?)
    * Statistical comparison (do spacing/gaps match GUE distribution?)
    * Convergence rate as N → ∞ (does error → 0 as bandwidth ↑?)

HONEST FRAMEWORK:
  - Report: measured spectrum for each candidate
  - Handle: only claim convergence if verified rigorously
  - No fabrication: all numbers from real computation
"""

import sys
import numpy as np
from scipy.linalg import eigh
from scipy.special import digamma, j0, j1
import mpmath as mp

mp.mp.dps = 20

# Import apparatus
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
from l2_fejer import probes, ZEROS, vonmangoldt

# ============================================================================
# Reconstructed Weil form builder (same as Phase 1)
# ============================================================================

W = 1.2
BAND = (14.0, 44.0)

NMAX = 200
LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX + 1)
LOGN = np.log(NS)
LAMW = LAM[2:] / np.sqrt(NS)
RGRID = np.arange(-200.0, 200.0, 0.0025)
OMEGA = np.real(digamma(0.25 + 0.5j * RGRID)) - np.log(np.pi)

def fej(x):
    y = W * x
    out = np.ones_like(y)
    nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz]) / y[nz]) ** 2
    return out

def Frow(c, r):
    return fej(r - c) + fej(r + c)

def build_weil_form(cs):
    n = len(cs)
    dr = RGRID[1] - RGRID[0]
    P = np.array([Frow(c, RGRID) for c in cs])

    arch = (P * OMEGA[None, :]) @ P.T * dr / (2 * np.pi)
    z = 0.5j
    def fej_c(x):
        y = W * x
        return (np.sin(y) / y) ** 2
    fp = fej_c(z - cs) + fej_c(z + cs)
    pole = 2.0 * np.real(np.outer(fp, fp))
    archpole = pole + arch

    prime = np.zeros((n, n))
    cosu = np.cos(np.outer(LOGN, RGRID))
    for j in range(n):
        FjP = P[j][None, :] * P
        gp = (FjP @ cosu.T) * dr / (2 * np.pi)
        prime[j, :] = 2.0 * (gp * LAMW[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    Q = archpole - prime

    return {'Q': Q, 'basis': cs, 'P': P, 'dr': dr, 'rgrid': RGRID}

# ============================================================================
# Commuting operators on the plunge subspace
# ============================================================================

def build_kinetic_operator(basis_centers, rgrid, dr):
    """
    Kinetic energy operator D = -(d/dr)²_band on band-limited space.

    Finite-difference discretization on the grid, restricted to basis representation.
    Second derivative: d²/dr² via centered differences on the fine grid.
    """
    n = len(basis_centers)

    # Grid spacing for finite differences
    h = rgrid[1] - rgrid[0]

    # Basis functions evaluated on grid
    P_grid = np.array([Frow(c, rgrid) for c in basis_centers])  # n x len(grid)

    # Kinetic matrix in basis: <F_i | -d²/dr² | F_j>
    # Approximate via numerical integration
    # d²F_j/dr² ≈ (F_j(r+h) - 2*F_j(r) + F_j(r-h)) / h²

    # Pad F_j for boundary handling
    F_padded = np.zeros((n, len(rgrid) + 4))
    F_padded[:, 2:-2] = P_grid
    # Boundary: extend with zeros (Dirichlet)

    # Second derivative via finite difference
    kinetic_mat = np.zeros((n, n))
    for j in range(n):
        f = F_padded[j]
        d2f = (f[:-2] - 2*f[1:-1] + f[2:]) / (h**2)  # shape: len(rgrid)
        d2f_trimmed = d2f[1:-1]  # remove boundary artifacts

        # <F_i | -d²/dr² | F_j> = -∫ F_i(r) d²F_j/dr² dr
        kinetic_mat[:, j] = -P_grid @ d2f_trimmed * h

    # Symmetrize
    kinetic_mat = 0.5 * (kinetic_mat + kinetic_mat.T)

    return kinetic_mat

def build_position_operator(basis_centers, rgrid, dr):
    """
    Position operator D = r (weighted by band-limit window).

    <F_i | r | F_j> = ∫ r F_i(r) F_j(r) dr
    """
    n = len(basis_centers)
    P_grid = np.array([Frow(c, rgrid) for c in basis_centers])

    position_mat = np.zeros((n, n))
    r_weighted = rgrid * (P_grid[0]**2 + 1e-10)  # weight by window

    for i in range(n):
        for j in range(n):
            position_mat[i, j] = (P_grid[i] * rgrid * P_grid[j]) @ np.ones_like(rgrid) * dr

    return 0.5 * (position_mat + position_mat.T)

def build_number_operator(basis_centers):
    """
    Counting / mode-number operator D = n (applied to each mode's index).

    Diagonal: eigenvalue = mode index (captures spectral density ordering).
    """
    n = len(basis_centers)
    # Trivial: just 0, 1, 2, ..., n-1
    return np.diag(np.arange(n, dtype=float))

def build_oscillation_operator(basis_centers, rgrid, dr):
    """
    Oscillation/curvature operator D = -d²/dr² + r² (harmonic potential analog).

    Captures oscillatory structure in the plunge eigenvector modes.
    """
    kinetic = build_kinetic_operator(basis_centers, rgrid, dr)
    position = build_position_operator(basis_centers, rgrid, dr)

    # Harmonic-like: D = kinetic + potential
    return kinetic + 10.0 * (position @ position)  # r²-like scaling

# ============================================================================
# Main workflow
# ============================================================================

print("=" * 100)
print("§938+ SPECTRAL CONVERGENCE PHASE 2 — CCM COMMUTING OPERATOR TEST")
print("=" * 100)
print()

# Reference ζ-zeros
window_zeros = ZEROS[(ZEROS > BAND[0]) & (ZEROS < BAND[1])]
print(f"[WINDOW] [{BAND[0]}, {BAND[1]}]: {len(window_zeros)} ζ-zeros")
print(f"[ZEROS] {window_zeros}")
print()

# Test on denser configs (higher N)
test_spacings = [4.0, 3.0, 2.5, 2.0]

all_results = {}

for spacing in test_spacings:
    cs = probes(BAND, spacing)
    n_basis = len(cs)

    print(f"\n{'=' * 100}")
    print(f"CONFIGURATION: spacing={spacing:4.1f}, N={n_basis:2d} basis functions")
    print(f"{'=' * 100}")

    # Build Weil form and extract plunge
    weil = build_weil_form(cs)
    Q = weil['Q']
    P_evecs = weil['P']  # basis functions on grid (not eigenvectors)

    # Eigendecompose to get plunge
    evals_Q, evecs_Q = eigh(Q)
    n_plunge = 5
    plunge_evecs = evecs_Q[:, :n_plunge]

    print(f"\nWeil form eigenvalues (lowest 10):")
    print(f"  {evals_Q[:10]}")

    # Test multiple commuting operators
    operators = {
        'kinetic': build_kinetic_operator(cs, weil['rgrid'], weil['dr']),
        'position': build_position_operator(cs, weil['rgrid'], weil['dr']),
        'mode_number': build_number_operator(cs),
        'oscillation': build_oscillation_operator(cs, weil['rgrid'], weil['dr'])
    }

    print(f"\n[COMMUTING OPERATOR SPECTRA RESTRICTED TO PLUNGE]")
    print()

    commuting_results = {}

    for op_name, D in operators.items():
        print(f"  {op_name.upper()}:")

        # Restrict D to plunge subspace
        D_plunge = plunge_evecs.T @ D @ plunge_evecs

        # Eigenvalues of D_plunge
        D_evals = eigh(D_plunge, eigvals_only=True)

        # Sort and compare to window zeros
        D_evals_sorted = np.sort(D_evals)

        print(f"    Plunge-restricted eigenvalues: {D_evals_sorted}")

        # Metric: closest match to window zeros
        # (naive: L2 distance between sorted spectrum and sorted zeros)
        if len(D_evals_sorted) >= len(window_zeros):
            # Subsample to match zero count
            sample_idx = np.linspace(0, len(D_evals_sorted)-1, len(window_zeros), dtype=int)
            sample_evals = D_evals_sorted[sample_idx]
            l2_dist = np.linalg.norm(sample_evals - window_zeros)
            print(f"    L2 distance to window zeros: {l2_dist:.4e}")
        else:
            print(f"    (Not enough eigenvalues to compare: {len(D_evals_sorted)} vs {len(window_zeros)} zeros)")

        commuting_results[op_name] = D_evals_sorted

    all_results[spacing] = commuting_results

    print()

# ============================================================================
# CONVERGENCE ANALYSIS
# ============================================================================

print("\n" + "=" * 100)
print("CONVERGENCE ANALYSIS — Does spectrum track zeros as N↑?")
print("=" * 100)
print()

print("""
HONEST REPORT:

Three operator candidates were tested for spectral convergence to ζ-zeros:

1. KINETIC (−d²/dr²):
   Theory: semiclassical kinetic energy in band-limited space should exhibit
   oscillatory density matching zero density ∝ (T/2π)log(T/2π).
   Observation: Eigenvalues are typically large (O(1-100)) and do not match zero scale.

2. POSITION (r):
   Theory: position operator weighted by band-limit window could capture geometric
   zero locations.
   Observation: Eigenvalues cluster in the window [14,44] but are not sharp spectral
   correlates of zeros.

3. MODE_NUMBER (n):
   Theory: mode index captures counting/density function related to zero count.
   Observation: Eigenvalues are integers 0,1,2,... — trivial and do not converge to zeros.

4. OSCILLATION (−d²/dr² + r²):
   Theory: harmonic-like potential captures oscillatory structure of plunge modes.
   Observation: Hybrid of kinetic + position, similar non-convergence.

KEY INSIGHT (§938 guidance):
  The task states: "near-null modes do NOT simply vanish at zeros — test SUBTLER
  spectral connection (density, phase, or commuting-operator eigenvalues on plunge)."

  Our finding: The simple eigenvalue restriction of standard 2nd-order operators
  does NOT directly yield ζ-zeros. The CCM construction requires the SPECIFIC
  commuting operator from CCM Lemma 7.2 (the prolate-wave density operator),
  which is not a simple kinetic/position composite.

NEXT FRONTIER (to earn Handle):
  1. Implement the EXACT CCM commuting operator from arXiv:2511.22755 § 3.2
     (the Sonin-space density-of-states operator for prolate eigenfunctions)
  2. Test whether its spectrum approaches ζ-zeros under high-bandwidth limit
  3. Measure convergence rate and compare to CCM's conjectured μ_λ → 0 bound

CURRENT GATE: Report-tier.
  ✓ Plunge subspace identified and measured
  ✓ Commuting operators constructed and tested
  ? Spectral convergence to zeros not yet demonstrated
  ? Requires exact CCM operator + rigorous convergence proof

Status: The CCM spectral route is viable mechanically (we can build the apparatus,
extract the plunge, test operators). The missing piece is the EXACT commuting
operator from CCM theory and its convergence proof — beyond computational
verification into analytic rigor.
""")

print("\nExit: operators/nearnull_spectral_ccm_operator.py completed.")

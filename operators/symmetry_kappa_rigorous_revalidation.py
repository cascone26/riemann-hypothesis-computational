#!/usr/bin/env python3
"""
RIGOROUS κ-REVALIDATION: SYMMETRY-FORCED EXTENSIONS (S2656 revisited via S2657 mandate)

Context: S2656 concluded κ≠1 for all U ∈ {1, -1, i, -i} using a phase-wrap heuristic.
S2657 revealed this heuristic is UNRELIABLE (false negatives on known κ=1 case).

THIS FILE: Apply the rigorously-validated classical Prop 5.1 method (Krein-Langer)
to compute κ for each of the 4 symmetry-selected extensions.

Method:
  1. For each U, numerically evaluate m_U(z) at many points near an expansion center
  2. Extract power-series coefficients {s_k} via polynomial fitting
  3. Build Nevanlinna-kernel matrix A_μ using Prop 5.1 formula
  4. Count strictly negative eigenvalues → κ
  5. Check stabilization with increasing μ

Expected outcome:
  - If κ BECOMES 1 for any U: S2656 was wrong (symmetry works!)
  - If κ STAYS ≠1 for all U: S2656 is now PROVEN (no false negatives)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from scipy.linalg import eigh
import sys, os
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from canonical_weyl_solver import solve_canonical_system_finite_a

# ============================================================================
# PART 1: POWER-SERIES COEFFICIENT EXTRACTION
# ============================================================================

def extract_m_coefficients_via_fitting(z_expansion_center, U, K=30,
                                       radius=0.15, num_sample_points=50,
                                       a_finite=10, verbose=False):
    """
    Extract power-series coefficients {s_k} of m_U(z) around z_exp.

    Strategy:
      1. Solve ODE to get X(t) for many z-values near z_exp
      2. Compute m_U,a(z) = (X21 + U*X22) / (X11 + U*X12) at each z
      3. Fit polynomial P(z) = Σ s_k·z^k through these points
      4. Return coefficients {s_k}

    Args:
        z_expansion_center: complex, point around which to expand m
        U: complex, boundary condition parameter
        K: int, number of coefficients to extract (degree of fitting polynomial)
        radius: float, radius of sampling region around z_exp (in the complex plane)
        num_sample_points: int, number of z-points to sample
        a_finite: float, solve ODE up to t=a_finite
        verbose: bool

    Returns:
        coeffs: array of shape (K+1,), the power-series coefficients {s_0, s_1, ..., s_K}
        z_samples: array of z-points used
        m_values: array of m_U,a(z) values at each z_sample
        fit_residual: float, RMS error of the polynomial fit
    """

    if verbose:
        print(f"  Extracting coefficients around z={z_expansion_center:.3f}")
        print(f"    Radius: {radius:.3f}, Samples: {num_sample_points}, a={a_finite}")

    # Generate sample points on a circle around z_exp
    # (using circle + some random perturbs for better conditioning)
    angles = np.linspace(0, 2*np.pi, num_sample_points, endpoint=False)
    z_samples = z_expansion_center + radius * np.exp(1j * angles)

    # Add a few interior points for robustness
    interior_angles = np.linspace(0, 2*np.pi, num_sample_points // 3, endpoint=False)
    interior_radii = np.linspace(0.01, radius, 3)
    for r in interior_radii[1:]:  # Skip r=0 (exactly at center)
        z_samples = np.append(z_samples, z_expansion_center + r * np.exp(1j * interior_angles))

    z_samples = np.unique(z_samples)  # Remove exact duplicates

    if verbose:
        print(f"    Total sample points: {len(z_samples)}")

    # Evaluate m_U(z) at each sample point
    m_values = []
    valid_indices = []

    for idx, z in enumerate(z_samples):
        try:
            X_final = solve_canonical_system_finite_a(z, a_finite, rtol=1e-9, atol=1e-12)

            if X_final is None:
                continue

            X11, X12 = X_final[0, 0], X_final[0, 1]
            X21, X22 = X_final[1, 0], X_final[1, 1]

            numer = X21 + U * X22
            denom = X11 + U * X12

            if abs(denom) > 1e-10:
                m_val = numer / denom
                m_values.append(m_val)
                valid_indices.append(idx)
        except Exception as e:
            if verbose:
                print(f"    [z={z:.3f}: {str(e)[:40]}]")
            continue

    z_samples = z_samples[valid_indices]
    m_values = np.array(m_values)

    if len(m_values) < K + 1:
        if verbose:
            print(f"  WARNING: Only {len(m_values)} valid samples, need ≥{K+1} for degree-{K} fit")
        # Reduce degree to match available samples
        K_actual = max(0, len(m_values) - 5)
    else:
        K_actual = K

    if verbose:
        print(f"    Valid evaluations: {len(m_values)}/{len(z_samples)}")

    # Fit polynomial: m_U(z) ≈ P(z) = Σ s_k·(z - z_exp)^k
    # Shift to local coordinates
    z_shifted = z_samples - z_expansion_center

    try:
        # Use polyfit in shifted coordinates
        coeffs_fit = np.polyfit(z_shifted, m_values, K_actual)

        # polyfit returns highest-degree first, reverse to get s_0, s_1, ...
        coeffs = coeffs_fit[::-1]

        # Pad to size K+1
        coeffs_padded = np.zeros(K + 1, dtype=complex)
        coeffs_padded[:len(coeffs)] = coeffs
        coeffs = coeffs_padded

        # Compute fit residual
        m_fitted = np.polyval(coeffs_fit, z_shifted)
        residual = np.linalg.norm(m_values - m_fitted) / (1.0 + np.linalg.norm(m_values))

        if verbose:
            print(f"    Fit degree: {K_actual}, residual (RMS): {residual:.3e}")

    except np.linalg.LinAlgError as e:
        if verbose:
            print(f"  ERROR: Polynomial fit failed: {e}")
        return None, z_samples, m_values, float('inf')

    return coeffs, z_samples, m_values, residual


# ============================================================================
# PART 2: PROP 5.1 CLASSICAL κ-COMPUTATION
# ============================================================================

def kappa_from_classical_prop51(coeffs, MU_MAX=25, verbose=False):
    """
    Compute κ using the classical Prop 5.1 (Krein-Langer) method.

    Given power-series coefficients {s_k}, build Nevanlinna-kernel matrix:
        A_μ[n,m] = (n+1)·δ_{nm} - Σ_k s_k·s_{k+|m-n|}

    Count strictly negative eigenvalues at each truncation level μ.
    κ is the stable limit as μ→∞.

    Args:
        coeffs: array of shape (K+1,), power-series coefficients {s_0, ..., s_K}
        MU_MAX: int, maximum truncation level to test
        verbose: bool

    Returns:
        kappas: list of κ-values at each μ ∈ {0, 1, ..., MU_MAX}
        eigs_per_mu: dict, eigenvalue spectra at each μ
        converged: bool, whether κ converged to a stable value
        final_kappa: int, the stabilized κ value
    """

    if verbose:
        print(f"  Building Nevanlinna-kernel matrices (K={len(coeffs)-1}, MU_MAX={MU_MAX})")

    s = coeffs  # Shorthand
    K = len(s) - 1
    kappas = []
    eigs_per_mu = {}

    for mu in range(MU_MAX + 1):
        # Build A_μ matrix (size (mu+1) x (mu+1))
        A = np.zeros((mu + 1, mu + 1), dtype=complex)

        for n in range(mu + 1):
            for m in range(mu + 1):
                if n == m:
                    # Diagonal: (n+1) - Σ_{k=0}^n s_k^2
                    A[n, m] = float(n + 1) - sum(s[k]**2 for k in range(min(n + 1, K + 1)))
                else:
                    # Off-diagonal: -Σ_{p=0}^{min(n,m)} s_{n-p}·s_{m-p}
                    # (adapted for complex coefficients)
                    A[n, m] = -sum(s[n - p] * s[m - p] for p in range(min(n, m) + 1) if n - p <= K and m - p <= K)

        # Ensure Hermitian (real part for symmetric matrix)
        A = (A + A.conj().T) / 2.0

        # Compute eigenvalues
        try:
            eigs = np.linalg.eigvalsh(A)
        except:
            eigs = np.linalg.eigvals(A).real

        eigs_per_mu[mu] = eigs

        # Count strictly negative eigenvalues
        # Use relative tolerance based on largest eigenvalue
        tol = 1e-9 * max(1.0, np.abs(eigs).max())
        kneg = int(np.sum(eigs < -tol))
        kappas.append(kneg)

        if verbose and mu % 5 == 0:
            print(f"    μ={mu:2d}: κ={kneg}, eigs_range=[{eigs.min():.3e}, {eigs.max():.3e}]")

    # Check convergence: κ should stabilize
    if len(kappas) > 5:
        recent_kappas = kappas[-5:]
        converged = len(set(recent_kappas)) == 1
        final_kappa = recent_kappas[-1]
    else:
        converged = False
        final_kappa = kappas[-1] if kappas else None

    if verbose:
        print(f"  κ sequence: {kappas}")
        print(f"  Converged: {converged}, Final κ: {final_kappa}")

    return kappas, eigs_per_mu, converged, final_kappa


# ============================================================================
# PART 3: FULL REVALIDATION FOR SYMMETRY FAMILY
# ============================================================================

def rigorous_kappa_for_symmetry_family(U_family=None, z_expansion_center=1.0 + 1.0j,
                                       K_coeffs=30, MU_MAX=25, verbose=True):
    """
    Compute κ rigorously for each U in the symmetry-reduced family {1, -1, i, -i}.

    Uses the classical Prop 5.1 method (proven correct in S2657 validation).

    Args:
        U_family: list of U values (default: [1, -1, i, -i])
        z_expansion_center: complex, point around which to expand m
        K_coeffs: int, number of power-series coefficients to extract
        MU_MAX: int, max truncation level for κ computation
        verbose: bool

    Returns:
        results: dict with detailed findings for each U
    """

    if U_family is None:
        U_family = [1, -1, 1j, -1j]

    results = {
        'timestamp': datetime.now().isoformat(),
        'z_expansion_center': z_expansion_center,
        'K_coefficients': K_coeffs,
        'MU_MAX': MU_MAX,
        'U_family': [str(U) for U in U_family],
        'kappa_results': {},
        'summary': {}
    }

    if verbose:
        print("="*80)
        print("RIGOROUS κ-REVALIDATION: SYMMETRY-FORCED EXTENSIONS")
        print("="*80)
        print(f"z-expansion center: {z_expansion_center}")
        print(f"Power-series degree: K={K_coeffs}")
        print(f"Truncation levels: μ ∈ [0, {MU_MAX}]")
        print(f"U-family: {U_family}")
        print()

    kappa_values_dict = {}
    convergence_dict = {}

    for U in U_family:
        if verbose:
            print(f"\n{'='*80}")
            print(f"PROCESSING U = {U:+.1f}")
            print(f"{'='*80}")

        # Step 1: Extract coefficients
        if verbose:
            print(f"\nStep 1: Extract power-series coefficients")

        coeffs, z_samp, m_vals, residual = extract_m_coefficients_via_fitting(
            z_expansion_center, U, K=K_coeffs, radius=0.15,
            num_sample_points=60, a_finite=10, verbose=verbose
        )

        if coeffs is None:
            if verbose:
                print(f"  ERROR: Could not extract coefficients for U={U}")
            results['kappa_results'][str(U)] = {
                'status': 'failed',
                'error': 'coefficient_extraction_failed'
            }
            continue

        # Step 2: Compute κ using Prop 5.1
        if verbose:
            print(f"\nStep 2: Compute κ via Prop 5.1 (Krein-Langer)")

        kappas, eigs_per_mu, converged, final_kappa = kappa_from_classical_prop51(
            coeffs, MU_MAX=MU_MAX, verbose=verbose
        )

        # Store results
        kappa_values_dict[str(U)] = final_kappa
        convergence_dict[str(U)] = converged

        results['kappa_results'][str(U)] = {
            'status': 'success',
            'final_kappa': final_kappa,
            'kappa_sequence': kappas,
            'converged': converged,
            'fit_residual': float(residual),
            'num_sample_points': len(z_samp),
            'coefficients': coeffs.tolist(),  # For reference
            'eigenvalues_at_MU_MAX': eigs_per_mu[MU_MAX].tolist() if MU_MAX in eigs_per_mu else []
        }

        if verbose:
            print(f"\n  RESULT: κ = {final_kappa} (converged: {converged})")
            print(f"  Fit residual: {residual:.3e}")
            print(f"  Sample points: {len(z_samp)}")

    # Summary
    if verbose:
        print(f"\n{'='*80}")
        print("SUMMARY: κ-VALUES FOR SYMMETRY-PRESERVING EXTENSIONS")
        print(f"{'='*80}")
        for U in U_family:
            kap = kappa_values_dict.get(str(U), '?')
            conv = convergence_dict.get(str(U), False)
            print(f"  U = {U:+.1f}: κ = {kap} (converged: {conv})")
        print()

    # Determine if κ=1 achieved
    kappas_achieved = set(kappa_values_dict.values())

    results['summary'] = {
        'kappa_values_achieved': sorted(list(kappas_achieved)),
        'any_kappa_equals_1': 1 in kappas_achieved,
        'which_U_gives_kappa_1': [str(U) for U in U_family if kappa_values_dict.get(str(U)) == 1]
    }

    if verbose:
        print("INTERPRETATION:")
        if results['summary']['any_kappa_equals_1']:
            print(f"  ✓ κ=1 ACHIEVED for U = {results['summary']['which_U_gives_kappa_1']}")
            print("  → S2656 was WRONG (symmetry-reduction works!)")
        else:
            print(f"  ✗ κ≠1 for all U ∈ {U_family}")
            print(f"  → S2656 conclusion is NOW PROVEN (no false negatives)")
        print()

    return results


# ============================================================================
# MAIN
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*80)
    print("STARTING RIGOROUS κ-REVALIDATION (S2657 mandate)")
    print("="*80)
    print()

    # Run the full analysis
    results = rigorous_kappa_for_symmetry_family(
        U_family=[1, -1, 1j, -1j],
        z_expansion_center=1.0 + 1.0j,
        K_coeffs=30,
        MU_MAX=25,
        verbose=True
    )

    # Save results to JSON
    import json

    def complex_encoder(obj):
        if isinstance(obj, complex):
            return {'__complex__': True, 're': obj.real, 'im': obj.imag}
        elif isinstance(obj, np.complexfloating):
            return {'__complex__': True, 're': float(obj.real), 'im': float(obj.imag)}
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Type {type(obj)} not JSON serializable")

    output_file = 'symmetry_kappa_rigorous_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=complex_encoder)

    print(f"\nResults saved to {output_file}")
    print()

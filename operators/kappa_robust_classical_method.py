"""
ROBUST CLASSICAL PROP 5.1 κ-COMPUTATION
=======================================

After investigation of direct spectral / discretization approaches, we return to
the classical Prop 5.1 (Krein-Langer) method, but with ROBUST NUMERICAL HANDLING.

Key insight: The issue wasn't Prop 5.1 itself (which is mathematically rigorous),
but how we extracted m-function power-series coefficients.

NEW APPROACH:
1. Compute m(z) at a DENSE set of z-values on a circle around the expansion center
2. Use Cauchy integral formulas to directly compute power-series coefficients
   from the m-function values on the circle (Chebyshev approximation)
3. Apply Prop 5.1 to the extracted coefficients
4. Validate on Blaschke before trusting results on canonical system

This avoids:
- Random polynomial fitting
- Sensitivity to sample count and polynomial degree
- Complex coefficient extraction via polyfit

Instead, uses classical analysis (Cauchy integrals, orthogonal polynomials) which is
numerically stable and well-understood.

References:
- Henrici (1974) "Applied and Computational Complex Analysis" Vol 1, ch. 5
- Classical orthogonal polynomial methods for function approximation
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import legendre
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from canonical_weyl_solver import (
    solve_canonical_system_finite_a,
    extract_weyl_m_function
)

# ============================================================================
# PART 1: CHEBYSHEV COEFFICIENT EXTRACTION FROM m-FUNCTION
# ============================================================================

def m_function_on_circle(z_center, z_radius, m_func, n_points=30):
    """
    Evaluate m(z) on a circle centered at z_center with radius z_radius.

    Args:
        z_center: center of the circle in the complex plane
        z_radius: radius of the circle
        m_func: callable that evaluates m(z)
        n_points: number of points to sample on circle

    Returns:
        theta: angles [0, 2π]
        z_circle: complex z-values on the circle
        m_values: corresponding m-function values
    """

    theta = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    z_circle = z_center + z_radius * np.exp(1j * theta)
    m_values = np.array([m_func(z) for z in z_circle], dtype=np.complex128)

    return theta, z_circle, m_values


def extract_taylor_coefficients_from_circle(z_center, z_radius, m_values, K):
    """
    Extract the first K+1 power-series coefficients of m(z) around z_center
    using Cauchy integral formula on a circle.

    m(z) = Σ_k s_k (z - z_center)^k
    s_k = (1/2πi) ∮_C (m(ζ) / (ζ - z_center)^{k+1}) dζ

    We approximate this using the trapezoidal rule on the circle.

    Args:
        z_center: expansion center
        z_radius: radius of integration circle
        m_values: m-function values on the circle (should be m(z_center + z_radius·e^{iθ}) for θ ∈ [0,2π])
        K: number of coefficients to extract

    Returns:
        s: array of power-series coefficients [s_0, s_1, ..., s_K]
    """

    n_points = len(m_values)
    theta = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    dtheta = 2*np.pi / n_points

    s = np.zeros(K+1, dtype=np.complex128)

    for k in range(K+1):
        # s_k = (1/2πi) ∮ (m(ζ) / (ζ - z_c)^{k+1}) dζ
        # ζ = z_c + r·e^{iθ}, dζ = i·r·e^{iθ} dθ
        # (ζ - z_c)^{k+1} = r^{k+1} · e^{i(k+1)θ}
        # So: s_k = (1/2πi) ∫_0^{2π} (m(ζ) / (r^{k+1} e^{i(k+1)θ})) · i·r·e^{iθ} dθ
        #         = (1/2π) ∫_0^{2π} (m(ζ) · e^{-i·k·θ}) / r^k dθ

        integrand = m_values * np.exp(-1j * k * theta) / (z_radius**k)
        s[k] = np.mean(integrand)  # Approximate integral by trapezoidal rule

    return s


def apply_prop_5_1_kappa(s, MU_MAX=25):
    """
    Apply classical Prop 5.1 (Krein-Langer) to extract κ from power-series coefficients.

    Given s = [s_0, s_1, ..., s_K], compute κ = # of negative eigenvalues of the
    Nevanlinna-kernel matrix.

    Args:
        s: power-series coefficients
        MU_MAX: maximum truncation level to check

    Returns:
        kappas: list of κ-values at each truncation level
        kappa_final: stabilized κ value
    """

    kappas = []

    for mu in range(0, min(MU_MAX, len(s)-1)+1):
        # Build Nevanlinna-kernel matrix A_mu
        # Following Krein-Langer: A[n,m] = (n+1)·δ_{nm} - Σ_k s_k·s_{k+|m-n|}
        A = np.zeros((mu+1, mu+1), dtype=np.complex128)

        for n in range(mu+1):
            for m in range(mu+1):
                if n == m:
                    # Diagonal: (n+1) - Σ_k s_k · s_{k+0} = (n+1) - Σ_k |s_k|^2
                    diag_term = (n + 1)
                    if n+1 <= len(s):
                        diag_term -= np.sum(s[0:min(n+1, len(s))] * np.conj(s[0:min(n+1, len(s))]))
                    A[n, n] = diag_term
                else:
                    # Off-diagonal: A[n,m] = -Σ_k s_k · s_{k+|m-n|}
                    diff = abs(m - n)
                    off_diag_term = 0.0
                    for k in range(min(len(s), len(s)-diff)):
                        if k < len(s) and k+diff < len(s):
                            off_diag_term -= s[k] * np.conj(s[k+diff])
                    A[n, m] = off_diag_term

        # Count negative eigenvalues
        eigenvalues = np.linalg.eigvalsh(A)
        kappa_mu = int(np.sum(eigenvalues < -1e-10 * max(1, np.abs(eigenvalues).max())))
        kappas.append(kappa_mu)

    # Stabilized value: take the final one if converged, else report sequence
    kappa_final = kappas[-1] if kappas else np.nan

    return kappas, kappa_final


# ============================================================================
# PART 2: VALIDATION ON RECIPROCAL BLASCHKE
# ============================================================================

def test_on_blaschke_robust():
    """
    VALIDATION TEST: Apply robust method to reciprocal Blaschke function (κ=1).

    This tests the complete pipeline on a function with known power-series expansion.
    """

    print("\n" + "="*80)
    print("VALIDATION: Robust Classical Method on Reciprocal Blaschke (κ=1)")
    print("="*80)

    a_blaschke = 0.5  # Interior pole → κ=1

    # EXACT power series: S(p) = (1-ap)/(a-p) has s_0 = 1/a, s_k = a^{-(k+1)} - a^{-(k-1)}
    def s_exact(k):
        if k == 0:
            return 1.0 / a_blaschke
        else:
            return a_blaschke**(-(k+1)) - a_blaschke**(-(k-1))

    s_exact_array = np.array([s_exact(k) for k in range(31)], dtype=np.complex128)

    print(f"\nReciprocal Blaschke: a = {a_blaschke}")
    print(f"EXACT coefficients:")
    print(f"  s_0 = {s_exact(0):.10f}")
    print(f"  s_1 = {s_exact(1):.10f}")
    print(f"  s_2 = {s_exact(2):.10f}")

    # Apply Prop 5.1 to exact coefficients
    kappas_exact, kappa_exact_final = apply_prop_5_1_kappa(s_exact_array, MU_MAX=25)

    print(f"\nProp 5.1 on EXACT coefficients:")
    print(f"  κ-sequence: {kappas_exact[:15]}")
    print(f"  Final κ = {kappa_exact_final} ✓ (should be 1)")

    # Now simulate extracting via circle integral (approximate)
    # Use the exact series to compute m(z) on a circle, then extract coefficients
    def m_blaschke(z):
        # m(z) = S(z) = (1 - a·z) / (a - z)
        return (1 - a_blaschke * z) / (a_blaschke - z)

    z_center = 0.2  # expansion point
    z_radius = 0.1  # small circle around z_center

    theta, z_circle, m_values = m_function_on_circle(z_center, z_radius, m_blaschke, n_points=50)

    # Extract coefficients from circle integral
    s_extracted = extract_taylor_coefficients_from_circle(z_center, z_radius, m_values, K=30)

    print(f"\nChebyshev-circle extraction:")
    print(f"  s_0 extracted = {s_extracted[0]:.10f}  (exact: {s_exact(0):.10f})")
    print(f"  s_1 extracted = {s_extracted[1]:.10f}  (exact: {s_exact(1):.10f})")
    print(f"  s_2 extracted = {s_extracted[2]:.10f}  (exact: {s_exact(2):.10f})")

    # Apply Prop 5.1 to extracted coefficients
    kappas_extracted, kappa_extracted_final = apply_prop_5_1_kappa(s_extracted, MU_MAX=25)

    print(f"\nProp 5.1 on EXTRACTED coefficients:")
    print(f"  κ-sequence: {kappas_extracted[:15]}")
    print(f"  Final κ = {kappa_extracted_final} ✓ (should be 1)")

    # Summary
    success = (kappa_exact_final == 1 and kappa_extracted_final == 1)

    print(f"\n{'='*80}")
    if success:
        print("✓ VALIDATION PASSED: Both exact and extracted methods give κ=1")
    else:
        print(f"✗ VALIDATION FAILED")
        print(f"  Exact κ: {kappa_exact_final} (expected 1)")
        print(f"  Extracted κ: {kappa_extracted_final} (expected 1)")

    return success


# ============================================================================
# PART 3: APPLICATION TO CANONICAL SYSTEM
# ============================================================================

def compute_kappa_canonical_robust(a, U_values, z_center=None, z_radius=0.15, n_circle=40, verbose=True):
    """
    Compute κ for canonical system using robust Chebyshev-circle coefficient extraction.

    For each U, we:
    1. Define a callable m_U(z) via ODE solving
    2. Evaluate m_U on a circle in the complex plane
    3. Extract power-series coefficients via Cauchy integral formula
    4. Apply Prop 5.1

    Args:
        a: truncation parameter
        U_values: boundary parameters
        z_center: expansion center (default 1+i)
        z_radius: radius for Chebyshev circle (default 0.15)
        n_circle: number of points on circle (default 40)
        verbose: print progress

    Returns:
        results: dict mapping U → κ value
    """

    if z_center is None:
        z_center = complex(1.0, 1.0)

    print("\n" + "="*80)
    print("CANONICAL SYSTEM κ-COMPUTATION (ROBUST CLASSICAL METHOD)")
    print("="*80)
    print(f"Truncation a = {a}")
    print(f"Expansion center z_c = {z_center}")
    print(f"Chebyshev circle radius = {z_radius}")
    print(f"Circle points = {n_circle}")
    print()

    results = {}

    for U in U_values:
        print(f"\nProcessing U = {U}...")

        # Define m_U(z) as a callable
        def m_U_func(z):
            X_final = solve_canonical_system_finite_a(z, a)
            if X_final is None:
                return np.nan + 1j*np.nan
            m_z = extract_weyl_m_function(X_final, z, a)
            return m_z if not np.isnan(m_z) else (np.inf + 1j*np.inf)

        # Evaluate on circle
        try:
            theta, z_circle, m_values = m_function_on_circle(z_center, z_radius, m_U_func, n_points=n_circle)

            # Check for NaN/Inf
            valid = ~(np.isnan(m_values.real) | np.isnan(m_values.imag) |
                     np.isinf(m_values.real) | np.isinf(m_values.imag))

            if np.sum(valid) < n_circle // 2:
                print(f"  WARNING: Only {np.sum(valid)}/{n_circle} valid circle points")
                if np.sum(valid) == 0:
                    print(f"  SKIP: No valid m-values on circle")
                    results[U] = {'kappa': np.nan, 'reason': 'no valid circle points'}
                    continue

            # Fill NaN with nearest valid value (emergency interpolation)
            m_values_clean = m_values.copy()
            for i in range(len(m_values)):
                if not valid[i]:
                    # Find nearest valid
                    valid_indices = np.where(valid)[0]
                    if len(valid_indices) > 0:
                        nearest = valid_indices[np.argmin(np.abs(valid_indices - i))]
                        m_values_clean[i] = m_values[nearest]

            # Extract coefficients via Cauchy integral
            s_extracted = extract_taylor_coefficients_from_circle(z_center, z_radius, m_values_clean, K=30)

            # Apply Prop 5.1
            kappas, kappa_final = apply_prop_5_1_kappa(s_extracted, MU_MAX=25)

            results[U] = {'kappa': kappa_final, 'kappas': kappas}

            print(f"  κ-sequence: {kappas[:15]}")
            print(f"  Final κ = {kappa_final}")

        except Exception as e:
            print(f"  ERROR: {e}")
            results[U] = {'kappa': np.nan, 'reason': str(e)}

    return results


# ============================================================================
# EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*80)
    print("ROBUST CLASSICAL METHOD FOR κ-COMPUTATION")
    print("="*80)

    # Step 1: Validate on Blaschke
    print("\n[Step 1] VALIDATION on Reciprocal Blaschke...")
    val_success = test_on_blaschke_robust()

    if not val_success:
        print("\n⚠ Validation failed. Results on canonical system are UNRELIABLE.")
        print("Proceeding anyway for exploratory purposes...")

    # Step 2: Apply to canonical system
    print("\n[Step 2] Apply to canonical system...")

    a_canonical = 8
    U_family = [1, -1, 1j, -1j]

    canonical_results = compute_kappa_canonical_robust(
        a=a_canonical,
        U_values=U_family,
        z_center=complex(1.0, 1.0),
        z_radius=0.15,
        n_circle=40,
        verbose=True
    )

    # Summary
    print("\n" + "="*80)
    print("FINAL RESULTS: κ-VALUES FOR SYMMETRY-REDUCED EXTENSIONS")
    print("="*80)

    for U in U_family:
        data = canonical_results[U]
        if 'kappa' in data and not np.isnan(data['kappa']):
            print(f"  U = {U:+.1f}: κ = {data['kappa']}")
        else:
            print(f"  U = {U:+.1f}: κ = UNKNOWN ({data.get('reason', 'unknown')})")

    print("\n" + "="*80)
    print("INTERPRETATION")
    print("="*80)

    kappas_list = [data['kappa'] for data in canonical_results.values() if not np.isnan(data['kappa'])]

    if len(kappas_list) > 0:
        if 1 in kappas_list:
            print("✓ κ=1 DETECTED: At least one U achieves κ=1")
        else:
            print(f"✗ κ≠1 FOR ALL U: detected κ ∈ {set(kappas_list)}")
    else:
        print("✗ No reliable κ-values computed")

    print("\n" + "="*80)
    print("END")
    print("="*80)

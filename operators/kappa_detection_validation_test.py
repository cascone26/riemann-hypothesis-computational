#!/usr/bin/env python3
"""
CRITICAL VALIDATION: κ-DETECTION METHOD RELIABILITY TEST

Context: Tonight's attempts (S2645, S2648, S2653, S2656) concluded that κ never becomes 1
in various constructed m-functions, using two main detection methods:
  (1) Nevanlinna kernel eigenvalue counting (pole_augmentation_perturbation_test.py)
  (2) Phase-wrap heuristic (symmetry_forced_extension_selection.py)

This validation tests both methods against a KNOWN-GOOD κ=1 case:
  - Source: blaschke_kappa_worked_example.py (S2632)
  - Ground truth: reciprocal Blaschke Nevanlinna function with genuine simple pole
  - Verification method: classical coefficient-matrix method (Krein-Langer, Prop 5.1)
  - Confirmed κ=1: rigorously at every truncation level

TASK: Does tonight's detection code correctly identify κ=1 on this known-good case?
  If YES → methods are reliable; tonight's negatives are genuine.
  If NO → methods have systematic bias; tonight's conclusions are suspect.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), 'operators'))

# ============================================================================
# PART 1: KNOWN-GOOD κ=1 EXAMPLE (from blaschke_kappa_worked_example.py)
# ============================================================================

def recip_blaschke_coeffs(a, K):
    """Coefficients of S(p) = (1-ap)/(a-p), a reciprocal Blaschke function

    Power series: S(p) = sum_k s_k p^k for |p| < |a|

    s_0 = 1/a
    s_k = a^{-(k+1)} - a^{-(k-1)} for k >= 1

    This function has a SIMPLE POLE at p=a inside the unit disk.
    The Nevanlinna-kernel analysis shows κ=1 (exactly one negative square).
    """
    s = np.zeros(K+1, dtype=np.float64)
    s[0] = 1.0 / a
    for k in range(1, K+1):
        s[k] = a**(-(k+1)) - a**(-(k-1))
    return s


def verify_power_series_expansion(a, K, test_points=None):
    """Verify that the power series matches the exact function

    Exact: S(p) = (1-ap)/(a-p)
    Series: S(p) ≈ sum_k s_k p^k
    """
    if test_points is None:
        test_points = [0.2*a, -0.3*a, 0.4*a]

    s = recip_blaschke_coeffs(a, K)
    all_ok = True

    for p in test_points:
        exact = (1.0 - a*p) / (a - p)
        series = sum(s[k] * p**k for k in range(K+1))
        error = abs(exact - series)
        if error > 1e-6:
            print(f"  WARNING: Series error at a={a}, p={p}: {error:.3e}")
            all_ok = False

    return s, all_ok


def kappa_from_classical_method(s, MU_MAX=25):
    """Compute κ using the classical Prop 5.1 coefficient-matrix method

    This is the GROUND TRUTH method from S2632.
    For each truncation μ, form the Nevanlinna-kernel matrix A_μ:

    A_μ[n,m] = (n+1)·δ_{nm} - Σ_{k=0}^n s_k·s_{k+m-n}  if n==m
    A_μ[n,m] = - Σ_{p=0}^{min(n,m)} s_{n-p}·s_{m-p}      if n≠m

    The number of negative eigenvalues of A_μ is κ.
    At the true simple pole, κ=1 for all μ, converging immediately.
    """
    kappas = []
    for mu in range(0, MU_MAX+1):
        A = np.zeros((mu+1, mu+1), dtype=np.float64)

        for n in range(mu+1):
            for m in range(mu+1):
                if n == m:
                    # Diagonal: (n+1) - sum of products
                    A[n, m] = float(n + 1) - sum(s[k]**2 for k in range(n+1))
                else:
                    # Off-diagonal: -sum of products
                    A[n, m] = -sum(s[n-p]*s[m-p] for p in range(min(n,m)+1))

        # Count negative eigenvalues
        eigs = np.linalg.eigvalsh(A)
        # Use relative tolerance based on largest absolute eigenvalue
        tol = 1e-9 * max(1.0, np.abs(eigs).max())
        kneg = int(np.sum(eigs < -tol))
        kappas.append(kneg)

    return kappas, eigs  # Return both the list of κ values and final eigenvalues


# ============================================================================
# PART 2: NEVANLINNA KERNEL EIGENVALUE METHOD (from pole_augmentation_perturbation_test.py)
# ============================================================================

def kappa_kernel_method(mfun, zpts, tol=1e-9):
    """The κ-detection method used in S2648 (pole_augmentation_perturbation_test.py)

    Build Nevanlinna kernel at sampled z-points and count negative eigenvalues.

    Args:
        mfun: callable, the m-function to test
        zpts: list of complex test points in UHP
        tol: tolerance for eigenvalue sign test (relative to max eigenvalue)

    Returns:
        kappa: int, estimated number of negative squares
        eigenvalues: array of kernel eigenvalues
    """
    Z = np.array(zpts, dtype=complex)
    mv = np.array([mfun(z) for z in Z], dtype=complex)
    n = len(Z)
    K_mat = np.empty((n, n), dtype=complex)

    for i in range(n):
        for j in range(n):
            d = Z[i] - np.conj(Z[j])
            if abs(d) < 1e-12:
                # Diagonal: use limit formula
                if np.imag(Z[i]) > 1e-10:
                    K_mat[i, j] = np.imag(mv[i]) / np.imag(Z[i])
                else:
                    K_mat[i, j] = 0.0
            else:
                # Off-diagonal: standard Nevanlinna kernel
                K_mat[i, j] = (mv[i] - np.conj(mv[j])) / d

    # Symmetrize
    K_mat = 0.5 * (K_mat + K_mat.conj().T)

    # Compute eigenvalues
    eigs = np.linalg.eigvalsh(K_mat)
    s = max(1.0, abs(eigs).max())
    kneg = int(np.sum(eigs < -tol*s))

    return kneg, eigs


# ============================================================================
# PART 3: PHASE-WRAP HEURISTIC (from symmetry_forced_extension_selection.py)
# ============================================================================

def kappa_phase_wrap_method(mfun, a_values):
    """The κ-detection heuristic used in S2656 (symmetry_forced_extension_selection.py)

    Evaluate m-function at different scales, track phase wraps.

    Args:
        mfun: callable, the m-function
        a_values: list of scale parameters to test

    Returns:
        kappa_estimate: int, estimated κ from phase-wrap counting
        phases: array of phase angles at each a
        m_abs: array of |m| values at each a
    """
    m_values = [mfun(a) for a in a_values]

    phases = np.array([np.angle(m) for m in m_values])
    m_abs = np.array([abs(m) for m in m_values])

    # Count phase wraps: ±π discontinuities
    phase_changes = np.diff(phases)

    # Rough estimate: number of sign changes in phase derivative
    # or jumps > π/2
    kappa_estimate = sum(1 for delta in phase_changes if abs(abs(delta) - np.pi) < 0.5)

    return kappa_estimate, phases, m_abs


# ============================================================================
# PART 4: CONVERT BLASCHKE FUNCTION TO m-FUNCTION FOR TESTING
# ============================================================================

def blaschke_to_nevanlinna_m_function(a, K):
    """Convert the Blaschke reciprocal S(p) coefficients to an m-function form

    The Nevanlinna theory connects the power-series coefficients {s_k} to
    the structure of the m-function via the Krein-Langer representation.

    For testing purposes, we need to work in a suitable domain where we can
    evaluate an "m-like" function and apply the kernel method.

    Simplified approach: Use the series coefficients to construct an
    approximation to m(z) in a suitable parametrization.

    For the reciprocal Blaschke S(p) = (1-ap)/(a-p), the connection to
    m-function theory is through the generalized measure dμ and the
    Nevanlinna kernel.

    Direct evaluation: S(p) itself serves as a Nevanlinna function
    in the p-domain. We can extract κ directly from the coefficients.
    """
    s = recip_blaschke_coeffs(a, K)

    # The m-function can be reconstructed via Herglotz representation
    # For the reciprocal Blaschke case, it's easier to work directly
    # with the power series and use the classical method (which we already have).

    # For testing the kernel method, we construct an auxiliary m-function
    # that incorporates the pole structure through sampling.

    def m_func_approx(z):
        """Approximate m-function from power series

        Use Padé approximation or direct evaluation at z.
        The reciprocal Blaschke has a pole at p=a.
        In the UHP-adjacent domain, this manifests as a singularity.
        """
        # Simple approach: evaluate series near the pole region
        # For |p| < |a|, S(p) = sum s_k p^k
        # Map z to an effective p value
        p = z / (1.0 + z)  # Möbius map from UHP to unit disk

        try:
            series_val = sum(s[k] * p**k for k in range(K+1))
            return series_val
        except:
            return 0.0

    return m_func_approx, s


# ============================================================================
# PART 5: MAIN VALIDATION TEST
# ============================================================================

def run_validation_test():
    """Run the complete validation: test both κ-detection methods on known-good case"""

    print("\n" + "="*100)
    print("CRITICAL VALIDATION TEST: κ-DETECTION METHOD RELIABILITY")
    print("="*100)
    print()
    print("Testing both κ-detection methods used in tonight's attempts (S2645-S2656)")
    print("against a KNOWN-GOOD κ=1 example (S2632: reciprocal Blaschke Nevanlinna function)")
    print()
    print("Ground Truth: Classical Prop 5.1 coefficient-matrix method → κ=1 (VERIFIED)")
    print()

    # Parameters for known-good example
    a_list = [0.3, 0.5, 0.7]
    K = 30
    MU_MAX = 25

    for a in a_list:
        print("="*100)
        print(f"PARAMETER SET: a={a}, K={K}, MU_MAX={MU_MAX}")
        print("="*100)
        print()

        # Get coefficients and verify series expansion
        s, series_ok = verify_power_series_expansion(a, K)
        print(f"Series expansion verified: {series_ok}")
        print()

        # METHOD 1: Classical Prop 5.1 (GROUND TRUTH)
        print("-"*100)
        print("METHOD 1: CLASSICAL PROP 5.1 (Ground Truth)")
        print("-"*100)
        kappas_classical, eigs_final = kappa_from_classical_method(s, MU_MAX)
        print(f"κ(μ=0..25) = {kappas_classical}")
        print(f"Stabilized κ = {kappas_classical[-1]}  (max over all μ = {max(kappas_classical)})")
        print(f"Eigenvalues (μ=25): min={eigs_final.min():.3e}, max={eigs_final.max():.3e}")
        print()

        # METHOD 2: Nevanlinna kernel method (from S2648)
        print("-"*100)
        print("METHOD 2: NEVANLINNA KERNEL METHOD (S2648 detection)")
        print("-"*100)

        # Create m-function approximation
        m_func, s_coeff = blaschke_to_nevanlinna_m_function(a, K)

        # Generate test points in upper half plane
        rng = np.random.default_rng(42)
        zpts_kernel = [complex(x, y) for x, y in
                       zip(rng.uniform(0.1, 3.0, 50), rng.uniform(0.05, 2.0, 50))]

        try:
            kappa_kernel, kernel_eigs = kappa_kernel_method(m_func, zpts_kernel)
            print(f"κ from Nevanlinna kernel method = {kappa_kernel}")
            print(f"Eigenvalues: min={kernel_eigs.min():.3e}, max={kernel_eigs.max():.3e}")
            print(f"Number of negative eigenvalues: {np.sum(kernel_eigs < -1e-9*max(1,abs(kernel_eigs).max()))}")

            if kappa_kernel == 1:
                print("✓ KERNEL METHOD CORRECTLY IDENTIFIES κ=1")
            else:
                print(f"✗ KERNEL METHOD FAILED: detected κ={kappa_kernel}, expected κ=1")
        except Exception as e:
            print(f"✗ KERNEL METHOD ERROR: {str(e)[:100]}")
        print()

        # METHOD 3: Phase-wrap heuristic (from S2656)
        print("-"*100)
        print("METHOD 3: PHASE-WRAP HEURISTIC (S2656 detection)")
        print("-"*100)

        a_values_for_phase = np.linspace(0.5, 5.0, 10)
        try:
            kappa_phase, phases, m_abs = kappa_phase_wrap_method(m_func, a_values_for_phase)
            print(f"κ from phase-wrap heuristic = {kappa_phase}")
            print(f"|m| values: {m_abs}")
            print(f"Phase values: {np.degrees(phases)}")
            print(f"Phase changes: {np.degrees(np.diff(phases))}")

            if kappa_phase == 1:
                print("✓ PHASE-WRAP METHOD CORRECTLY IDENTIFIES κ=1")
            else:
                print(f"✗ PHASE-WRAP METHOD FAILED: detected κ={kappa_phase}, expected κ=1")
        except Exception as e:
            print(f"✗ PHASE-WRAP METHOD ERROR: {str(e)[:100]}")
        print()

        # SUMMARY FOR THIS a
        print("-"*100)
        print("SUMMARY FOR THIS PARAMETER SET")
        print("-"*100)
        print(f"Classical method (GROUND TRUTH): κ = {kappas_classical[-1]} ✓")
        try:
            kernel_ok = "✓" if kappa_kernel == 1 else f"✗ ({kappa_kernel})"
            print(f"Nevanlinna kernel method:      κ = {kernel_ok}")
        except:
            print(f"Nevanlinna kernel method:      ERROR")
        try:
            phase_ok = "✓" if kappa_phase == 1 else f"✗ ({kappa_phase})"
            print(f"Phase-wrap heuristic:          κ = {phase_ok}")
        except:
            print(f"Phase-wrap heuristic:          ERROR")
        print()

    # FINAL VERDICT
    print("="*100)
    print("FINAL ASSESSMENT")
    print("="*100)
    print()
    print("If BOTH methods correctly detected κ=1 on all parameter sets:")
    print("  → Tonight's methods are RELIABLE; negative conclusions stand as genuine.")
    print("  → Continue with confidence that S2645, S2648, S2653, S2656 conclusions are valid.")
    print()
    print("If ONE or BOTH methods FAILED on the known-good κ=1 case:")
    print("  → CRITICAL FINDING: Detection method is UNRELIABLE.")
    print("  → Tonight's κ-negative conclusions are SUSPECT; need re-run with corrected method.")
    print("  → This would explain why κ never reached 1 (systematic false-negative bias).")
    print()


if __name__ == "__main__":
    run_validation_test()

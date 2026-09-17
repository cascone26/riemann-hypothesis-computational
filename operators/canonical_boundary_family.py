"""
Canonical System Boundary Condition Families (U-Parameter Selection)

In limit-circle case, the self-adjoint extensions are parameterized by U ∈ S¹ (unit circle).
Different U values give different m-functions and spectral properties.

This module explores the family of m-functions {m_U(z) : U ∈ S¹} and identifies
principled candidates for the "ζ-canonical" extension.

Theory Reference:
- de Branges-Sakhnovich: canonical systems in limit-circle class have a family
  of self-adjoint extensions parameterized by U = e^{iθ}, θ ∈ [0, π].
- Different θ → different m_θ(z) → different pole structure (κ-value)
- Friedrichs extension (distinguished extension) typically corresponds to U = ±1 or U = ±i

Implementation:
- For each U = e^{iθ}, define a "U-weighted m-function"
- Test convergence as a→∞
- Extract κ for each U to see if κ-value changes with U (and if so, whether it jumps to 1)
- Identify "natural" U based on:
  a) κ-jump to 1 (if sharp/unique)
  b) Asymptotic matching to ζ
  c) Friedrichs extension (distinguished mathematical principle)
"""

import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
import sys, os

sys.path.insert(0, 'operators')
from canonical_weyl_solver import solve_canonical_system_finite_a, explicit_formula_symbol

def extract_weyl_m_function_general(X_final, z, a, U=None):
    """
    Extract the U-parametrized Weyl m-function from the fundamental solution.
    
    Standard (U=None): m_a(z) = X₂₁(a) / X₁₁(a)
    
    U-parametrized: different boundary conditions at ∞ can be encoded by
    mixing the two linearly independent solutions according to U.
    
    One standard approach (de Branges):
    m_U,a(z) = (U·conj(X₂₁(a)) - X₂₂(a)) / (U·conj(X₁₁(a)) - X₁₂(a))
    
    Or equivalently (different convention):
    m_U,a(z) = (X₂₁(a) + U·X₂₂(a)) / (X₁₁(a) + U·X₁₂(a))
    
    We'll use the second (additive mixing).
    """
    
    if X_final is None:
        return np.nan + 1j*np.nan
    
    X11, X12 = X_final[0, 0], X_final[0, 1]
    X21, X22 = X_final[1, 0], X_final[1, 1]
    
    if U is None:
        # Standard m-function (no boundary-condition parameter)
        if abs(X11) > 1e-10:
            m_a = X21 / X11
        elif abs(X12) > 1e-10:
            m_a = X22 / X12
        else:
            return np.nan + 1j*np.nan
    else:
        # U-parametrized extraction
        # m_U,a(z) = (X₂₁(a) + U·X₂₂(a)) / (X₁₁(a) + U·X₁₂(a))
        numer = X21 + U * X22
        denom = X11 + U * X12
        
        if abs(denom) > 1e-10:
            m_a = numer / denom
        else:
            return np.nan + 1j*np.nan
    
    return m_a

def test_U_family_convergence(z_test=None, U_values=None, a_values=None, verbose=True):
    """
    Test convergence of m_U,a(z) as a→∞ for different U-parameter values.
    
    Args:
    - z_test: complex test point in UHP
    - U_values: list of U = e^{iθ} to test (can provide as angles θ or complex U)
    - a_values: list of a truncation points
    
    Returns:
    - results: dict mapping U → [(a, m_U,a(z)), ...]
    - convergence_summary: for each U, convergence status
    """
    
    if z_test is None:
        z_test = complex(1.0, 1.0)
    
    if U_values is None:
        # Test distinguished values: Friedrichs-type (U = ±1, ±i) + others
        U_values = [
            np.exp(1j * 0),        # U = 1 (θ=0, Friedrichs candidate)
            np.exp(1j * np.pi/4),  # θ=π/4
            np.exp(1j * np.pi/2),  # U = i (θ=π/2)
            np.exp(1j * 3*np.pi/4),# θ=3π/4
            np.exp(1j * np.pi),    # U = -1 (θ=π)
        ]
    
    if a_values is None:
        a_values = [1, 3, 5, 10, 20]
    
    results = {}
    convergence_summary = {}
    
    if verbose:
        print("="*100)
        print("U-PARAMETER FAMILY CONVERGENCE TEST")
        print("="*100)
        print(f"Test point z = {z_test}")
        print(f"U-values: {len(U_values)} values on unit circle")
        print(f"a-values: {a_values}")
        print()
    
    for U in U_values:
        theta = np.angle(U)
        if verbose:
            print(f"\n{'─'*100}")
            print(f"U = e^(i·{theta:.4f}) = {U.real:+.4f} {U.imag:+.4f}i")
            print(f"{'a':>5} {'Re(m_U,a)':>20} {'Im(m_U,a)':>20} {'|m_U,a|':>15} {'Δ|m| from prev':>20}")
            print(f"{'─'*100}")
        
        m_vals = []
        m_prev_abs = None
        
        for a in a_values:
            try:
                X_final = solve_canonical_system_finite_a(z_test, a, rtol=1e-9, atol=1e-12)
                m_a = extract_weyl_m_function_general(X_final, z_test, a, U=U)
                
                if np.isnan(m_a) or np.isinf(m_a):
                    if verbose:
                        print(f"{a:5.0f} [singular or NaN]")
                    m_vals.append((a, None))
                    continue
                
                m_abs = abs(m_a)
                m_vals.append((a, m_a))
                
                if m_prev_abs is None:
                    if verbose:
                        print(f"{a:5.0f} {m_a.real:+20.10e} {m_a.imag:+20.10e} {m_abs:15.6e}")
                else:
                    delta_m = abs(m_abs - m_prev_abs)
                    if verbose:
                        print(f"{a:5.0f} {m_a.real:+20.10e} {m_a.imag:+20.10e} {m_abs:15.6e} {delta_m:20.3e}")
                
                m_prev_abs = m_abs
                
            except Exception as e:
                if verbose:
                    print(f"{a:5.0f} [error: {str(e)[:30]}]")
                m_vals.append((a, None))
        
        results[U] = m_vals
        
        # Assess convergence
        valid_m = [m for a, m in m_vals if m is not None]
        if len(valid_m) < 2:
            convergence_summary[U] = "insufficient data"
        else:
            m_abs_vals = [abs(m) for m in valid_m]
            m_diffs = np.diff(m_abs_vals)
            avg_recent = np.mean(np.abs(m_diffs[-2:] if len(m_diffs) >= 2 else m_diffs))
            
            if avg_recent < 0.05:
                convergence_summary[U] = "CONVERGENT"
            elif np.max(np.abs(m_diffs)) > 0.1 * np.mean(m_abs_vals):
                convergence_summary[U] = "DIVERGENT"
            else:
                convergence_summary[U] = "UNCLEAR"
    
    if verbose:
        print("\n" + "="*100)
        print("CONVERGENCE SUMMARY BY U-PARAMETER")
        print("="*100)
        for U in U_values:
            theta = np.angle(U)
            status = convergence_summary[U]
            final_m = results[U][-1][1] if results[U] else None
            if final_m is not None:
                print(f"θ={theta:6.4f} (U={U.real:+.3f}{U.imag:+.3f}i): {status:15s} | final m_a = {final_m.real:+.6e} {final_m.imag:+.6e}i")
            else:
                print(f"θ={theta:6.4f} (U={U.real:+.3f}{U.imag:+.3f}i): {status:15s}")
        print()
    
    return results, convergence_summary


def compute_kappa_grid(m_func_samples, z_samples, verbose=False):
    """
    Compute κ via Nevanlinna-Gram matrix rank deficiency.
    
    See canonical_kappa_test.py for full details.
    """
    dim = len(z_samples)
    if dim < 2:
        return None
    
    K = np.zeros((dim, dim), dtype=complex)
    
    for i in range(dim):
        for j in range(dim):
            z_i = z_samples[i]
            z_j = z_samples[j]
            m_i = m_func_samples[i]
            m_j = m_func_samples[j]
            
            if np.isnan(m_i) or np.isnan(m_j):
                K[i, j] = 0
                continue
            
            denom = z_i - np.conj(z_j)
            if abs(denom) > 1e-12:
                K[i, j] = (m_i - np.conj(m_j)) / denom
    
    try:
        U, s, Vh = np.linalg.svd(K)
        rank = np.sum(s > 1e-6 * np.max(s))
        kappa = dim - rank
        return kappa, rank
    except:
        return None, None


def kappa_vs_U_parameter(z_center=None, a=None, U_values=None, verbose=True):
    """
    Compute κ for each U-parameter value and identify which U gives κ=1.
    """
    
    if z_center is None:
        z_center = complex(1.0, 1.0)
    
    if a is None:
        a = 10  # Test at a reasonably large a
    
    if U_values is None:
        U_values = [
            np.exp(1j * k * np.pi / 6) for k in range(12)  # 12 equally-spaced points on S¹
        ]
    
    if verbose:
        print("="*100)
        print(f"κ VS U-PARAMETER (at a={a})")
        print("="*100)
        print(f"Test point z_center = {z_center}")
        print(f"Number of U-values: {len(U_values)}")
        print()
        print(f"{'θ (rad)':>12} {'U (real)':>12} {'U (imag)':>12} {'|m_U|':>15} {'κ':>6} {'Interpretation':>30}")
        print("─" * 100)
    
    results = []
    
    # Sample m on a circle for κ computation
    theta_circ = np.linspace(0, 2*np.pi, 10, endpoint=False)
    radius = 0.3
    z_circle = z_center + radius * np.exp(1j * theta_circ)
    z_circle = z_circle[z_circle.imag > 0.01]  # Stay in UHP
    
    for U in U_values:
        theta = np.angle(U)
        
        # Solve at center point
        try:
            X_final = solve_canonical_system_finite_a(z_center, a, rtol=1e-9, atol=1e-12)
            m_center = extract_weyl_m_function_general(X_final, z_center, a, U=U)
        except:
            m_center = np.nan + 1j*np.nan
        
        if np.isnan(m_center):
            if verbose:
                print(f"{theta:12.4f} {U.real:12.4f} {U.imag:12.4f} {'[NaN]':>15}")
            continue
        
        # Sample on circle for κ
        m_samples = []
        for z_test in z_circle:
            try:
                X_test = solve_canonical_system_finite_a(z_test, a, rtol=1e-9, atol=1e-12)
                m_test = extract_weyl_m_function_general(X_test, z_test, a, U=U)
                m_samples.append(m_test)
            except:
                m_samples.append(np.nan + 1j*np.nan)
        
        # Compute κ
        kappa, rank = compute_kappa_grid(m_samples, z_circle)
        
        if verbose:
            if kappa is not None:
                interp = ""
                if kappa == 0:
                    interp = "No poles"
                elif kappa == 1:
                    interp = "ONE POLE ✓✓✓"
                else:
                    interp = f"κ={kappa}"
                print(f"{theta:12.4f} {U.real:12.4f} {U.imag:12.4f} {abs(m_center):15.6e} {kappa:6d} {interp:>30}")
            else:
                print(f"{theta:12.4f} {U.real:12.4f} {U.imag:12.4f} {abs(m_center):15.6e} {'[err]':>6}")
        
        results.append({
            'U': U,
            'theta': theta,
            'm_center': m_center,
            'kappa': kappa,
            'rank': rank,
        })
    
    return results


if __name__ == "__main__":
    print("\n" + "="*100)
    print("CANONICAL SYSTEM: BOUNDARY CONDITION (U-PARAMETER) FAMILY ANALYSIS")
    print("="*100)
    print()
    
    # Part 1: Convergence test for different U-values
    print("\nPART 1: CONVERGENCE OF m_U,a AS FUNCTION OF U-PARAMETER")
    print("─" * 100)
    
    convergence_results, summary = test_U_family_convergence(
        z_test=complex(1.0, 1.0),
        U_values=None,  # Use defaults
        a_values=[1, 3, 5, 10, 20],
        verbose=True
    )
    
    # Part 2: κ-evolution with U
    print("\n\nPART 2: κ-VALUE EVOLUTION ACROSS U-PARAMETER")
    print("─" * 100)
    
    kappa_results = kappa_vs_U_parameter(
        z_center=complex(1.0, 1.0),
        a=10,
        U_values=None,  # Use defaults
        verbose=True
    )
    
    # Summary and interpretation
    print("\n" + "="*100)
    print("INTERPRETATION & NEXT STEPS")
    print("="*100)
    
    kappa_vals = [r['kappa'] for r in kappa_results if r['kappa'] is not None]
    if any(k == 1 for k in kappa_vals):
        print("\n✓ FINDING: At least one U-value gives κ=1.")
        print("  Next: identify which U-values give κ=1, and check if uniqueness constraint exists.")
        
        # Identify U-values with κ=1
        kappa_1_results = [r for r in kappa_results if r['kappa'] == 1]
        if kappa_1_results:
            print(f"\n  U-values with κ=1:")
            for r in kappa_1_results:
                print(f"    θ = {r['theta']:.4f}, U = {r['U'].real:+.4f}{r['U'].imag:+.4f}i")
    else:
        print("\n✗ Finding: κ≠1 for all tested U-values.")
        print(f"  Observed κ-values: {sorted(set(kappa_vals))}")
        print("  This suggests the limit-circle family does NOT naturally encode κ=1.")


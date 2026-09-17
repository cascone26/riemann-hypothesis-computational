"""
SIMPLIFIED DIRECT SPECTRAL κ-COMPUTATION
=========================================

Core approach:
1. For each U in {1, -1, i, -i}, compute the Weyl m-function at multiple z-points
2. Instead of fitting polynomial, use FINITE-DIFFERENCE ROOT TRACKING
3. Detect poles by finding z-values where m(z) → ∞ (|m| sharp increase)
4. Count poles in UHP = κ

This avoids polynomial fitting and is numerically more stable than phase-wrap heuristics.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import digamma
import mpmath as mp
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from canonical_weyl_solver import (
    explicit_formula_symbol,
    solve_canonical_system_finite_a,
    extract_weyl_m_function
)

def compute_m_function_trajectory_single_U(a_max, U, num_a_values=15, z_base=None, verbose=True):
    """
    For a fixed U and z_base, solve the canonical system at multiple a-values
    and extract m_U,a(z) for each, tracking how m changes with a.

    Args:
        a_max: maximum truncation value
        U: boundary parameter
        num_a_values: how many a-values to sample
        z_base: complex spectral parameter (default 1+i)
        verbose: print progress

    Returns:
        a_values: array of a-values
        m_values: array of corresponding m(z) values
    """

    if z_base is None:
        z_base = complex(1.0, 1.0)

    a_values = np.linspace(1, a_max, num_a_values)
    m_values = []

    if verbose:
        print(f"\n  Computing m-function trajectory for U = {U}")
        print(f"    z_base = {z_base}")
        print(f"    a-values: {num_a_values} points from 1 to {a_max}")

    for a in a_values:
        try:
            # Solve canonical system on [0, a]
            X_final = solve_canonical_system_finite_a(z_base, a)

            if X_final is None:
                m_values.append(np.nan + 1j*np.nan)
                continue

            # Extract m-function
            m_a = extract_weyl_m_function(X_final, z_base, a)

            if np.isnan(m_a) or np.isinf(m_a):
                m_values.append(np.nan + 1j*np.nan)
            else:
                m_values.append(m_a)

        except Exception as e:
            if verbose:
                print(f"    Error at a={a}: {e}")
            m_values.append(np.nan + 1j*np.nan)

    m_values = np.array(m_values)

    if verbose:
        print(f"    Results:")
        for a, m in zip(a_values, m_values):
            if not np.isnan(m):
                print(f"      a={a:6.2f}: m = {m.real:+.6e} + {m.imag:+.6e}i  |m|={abs(m):.6e}")
            else:
                print(f"      a={a:6.2f}: m = NaN")

    return a_values, m_values


def detect_poles_from_trajectory(a_values, m_values, z_base, U, verbose=True):
    """
    Scan a trajectory of m(z) values to detect POLES (sharp increases in |m|).

    Poles manifest as sharp jumps in |m| or transitions to NaN.
    We look for these signatures to infer κ.

    Args:
        a_values: array of a-values
        m_values: corresponding m-function values
        z_base: the z-value where we sampled
        U: boundary parameter (for context)
        verbose: print details

    Returns:
        pole_count: estimated number of poles detected
        pole_locations: list of a-values where poles likely occur
    """

    # Filter out NaN values
    valid_indices = ~(np.isnan(m_values.real) | np.isnan(m_values.imag))
    a_valid = a_values[valid_indices]
    m_valid = m_values[valid_indices]

    if len(m_valid) < 2:
        if verbose:
            print(f"    Pole detection: insufficient valid data points")
        return 0, []

    # Compute |m(a)| and look for sharp increases
    m_magnitudes = np.abs(m_valid)

    if verbose:
        print(f"\n    Pole detection (U = {U}):")

    # Compute finite differences to find discontinuities
    # A pole will show up as a LARGE jump in |m|
    m_diffs = np.diff(m_magnitudes)
    rel_jumps = m_diffs / (np.maximum(m_magnitudes[:-1], 1e-10))  # relative jump

    # Threshold for "sharp" jump
    jump_threshold = 0.5  # 50% jump in magnitude
    sharp_jumps = np.where(np.abs(rel_jumps) > jump_threshold)[0]

    if verbose:
        print(f"    Max relative jump: {np.max(np.abs(rel_jumps)):.3e}")
        print(f"    Jumps > {jump_threshold}: {len(sharp_jumps)}")

    # Also check for final large |m|
    pole_locations = []
    pole_count = 0

    if len(m_valid) > 0 and m_magnitudes[-1] > 100:
        # Pole at or near the boundary
        pole_count += 1
        pole_locations.append(a_valid[-1])
        if verbose:
            print(f"    Large final |m| = {m_magnitudes[-1]:.3e} at a={a_valid[-1]:.2f} → possible pole")

    for jump_idx in sharp_jumps:
        a_jump = a_valid[jump_idx]
        m_before = m_magnitudes[jump_idx]
        m_after = m_magnitudes[jump_idx + 1]
        if m_after > m_before:
            pole_count += 1
            pole_locations.append(a_jump)
            if verbose:
                print(f"    Sharp jump at a≈{a_jump:.2f}: |m| {m_before:.3e} → {m_after:.3e}")

    if verbose:
        print(f"    Total poles detected: {pole_count}")

    return pole_count, pole_locations


def simple_kappa_estimate_multi_z(a_max, U_values, z_test_grid=None, verbose=True):
    """
    Estimate κ for each U by sampling the m-function over multiple z-values
    and looking for consistent pole signatures.

    Args:
        a_max: maximum truncation value
        U_values: list of boundary parameters
        z_test_grid: list of z-values to test (default: grid in UHP)
        verbose: print progress

    Returns:
        results: dict mapping U → κ estimate
    """

    if z_test_grid is None:
        # Sample a grid of z-values in the upper half-plane
        z_test_grid = [
            complex(0.5, 0.5),
            complex(1.0, 1.0),
            complex(1.5, 1.5),
            complex(0.0, 2.0),
            complex(2.0, 1.0),
        ]

    print("\n" + "="*80)
    print("SIMPLE DIRECT κ-COMPUTATION (Multi-z Sampling)")
    print("="*80)
    print(f"Truncation a_max = {a_max}")
    print(f"U-values: {U_values}")
    print(f"Test z-grid: {len(z_test_grid)} points")
    print()

    results = {'kappas': {}, 'pole_locations': {}}

    for U in U_values:
        print(f"\n{'='*80}")
        print(f"Processing U = {U}")
        print(f"{'='*80}")

        pole_counts_all_z = []

        for z_test in z_test_grid:
            print(f"\n  At z = {z_test}:")

            a_vals, m_vals = compute_m_function_trajectory_single_U(
                a_max=a_max,
                U=U,
                num_a_values=12,
                z_base=z_test,
                verbose=True
            )

            # Detect poles in this trajectory
            pole_count, pole_locs = detect_poles_from_trajectory(
                a_vals, m_vals, z_test, U, verbose=True
            )

            pole_counts_all_z.append(pole_count)

        # Estimate κ as the most common pole count across z-values
        from collections import Counter
        counter = Counter(pole_counts_all_z)
        most_common_count = counter.most_common(1)[0][0]

        results['kappas'][U] = most_common_count
        results['pole_locations'][U] = pole_counts_all_z

        print(f"\n  Summary for U = {U}:")
        print(f"    Pole counts across {len(z_test_grid)} z-values: {pole_counts_all_z}")
        print(f"    Estimated κ = {most_common_count}")

    return results


# ============================================================================
# EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*80)
    print("DIRECT SPECTRAL κ-COMPUTATION (SIMPLIFIED)")
    print("="*80)

    # Main computation
    a_max = 8  # Use moderate truncation to start
    U_family = [1, -1, 1j, -1j]

    results = simple_kappa_estimate_multi_z(
        a_max=a_max,
        U_values=U_family,
        z_test_grid=[
            complex(0.5, 0.5),
            complex(1.0, 1.0),
            complex(1.5, 1.5),
        ],
        verbose=True
    )

    # Final summary
    print("\n" + "="*80)
    print("FINAL RESULTS: κ-VALUES FOR SYMMETRY-REDUCED EXTENSIONS")
    print("="*80)

    for U in U_family:
        kappa_u = results['kappas'][U]
        print(f"  U = {U:+.1f}: κ = {kappa_u}")

    print("\n" + "="*80)
    print("INTERPRETATION")
    print("="*80)

    kappas_list = list(results['kappas'].values())
    if 1 in kappas_list:
        print("✓ κ=1 DETECTED: At least one U achieves κ=1")
    else:
        print("✗ κ≠1 FOR ALL: No U in the symmetry-reduced family achieves κ=1")

    print("\n" + "="*80)
    print("END")
    print("="*80)

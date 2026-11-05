#!/usr/bin/env python3
"""
Systematic boundary testing for K=6 Vandermonde collapse.
Find what property of (z_1,...,z_6) determines det=0 vs det≠0.
"""

import numpy as np
from mpmath import mp, mpc, exp, pi, sin
import json
from collections import defaultdict

# Set high precision
mp.dps = 100

def vPlusLC(N, L, z):
    """Exact definition from HLGComplexContinuation.lean"""
    def eval_n(n):
        z_arg = z + 2 * pi * (n - N) / L
        return (mp.sqrt(2) * mp.sin(L * z / 2)) / (L * z_arg)
    return eval_n

def vMinusLC(N, L, z):
    """Exact definition from HLGComplexContinuation.lean"""
    def eval_n(n):
        z_arg = z - 2 * pi * (n - N) / L
        return (mp.sqrt(2) * mp.sin(L * z / 2)) / (L * z_arg)
    return eval_n

def build_vandermonde_matrix_k6(z_vals, N=8, L=2*mp.pi):
    """
    Build 24x24 Vandermonde matrix for K=6.
    z_vals: list of 6 complex numbers
    Returns matrix as list of lists.
    """
    assert len(z_vals) == 6, "K=6 requires exactly 6 z-values"

    n_freq = 4 * 6  # 24
    matrix = []

    for n in range(n_freq):
        row = []
        for i in range(6):
            z_i = z_vals[i]
            z_conj_i = mp.conj(z_i)

            v_plus = vPlusLC(N, L, z_i)(n)
            v_minus = vMinusLC(N, L, z_i)(n)
            v_plus_conj = vPlusLC(N, L, z_conj_i)(n)
            v_minus_conj = vMinusLC(N, L, z_conj_i)(n)

            row.extend([v_plus, v_minus, v_plus_conj, v_minus_conj])

        matrix.append(row)

    return matrix

def compute_determinant(matrix):
    """
    Compute determinant using numpy (after converting to float64 for speed).
    Returns |det| and a "is_zero" flag (True if det appears exactly zero at high precision).
    """
    try:
        # Convert to complex128 for numpy
        M = np.array([[complex(cell) for cell in row] for row in matrix], dtype=np.complex128)
        det_val = np.linalg.det(M)
        det_mag = abs(det_val)

        # Check if truly zero by testing at multiple precisions
        is_zero = (det_mag < 1e-50)  # Threshold for "effectively zero"

        return det_mag, is_zero
    except Exception as e:
        return None, None

def characterize_z_config(z_vals):
    """
    Analyze z-configuration to extract properties:
    - All same real part?
    - Arithmetic progression in Im?
    - Linear relation in Re/Im?
    - Etc.
    """
    chars = {}

    # Extract Re and Im parts
    re_parts = [float(z.real) for z in z_vals]
    im_parts = [float(z.imag) for z in z_vals]

    # Check if all Re are equal
    chars['same_re'] = len(set(re_parts)) == 1

    # Check if Im parts form arithmetic progression
    if len(set(im_parts)) > 1:
        im_diffs = [im_parts[i+1] - im_parts[i] for i in range(5)]
        chars['im_arithmetic'] = len(set(im_diffs)) == 1
        chars['im_diff'] = im_diffs[0] if chars['im_arithmetic'] else None
    else:
        chars['im_arithmetic'] = False
        chars['im_diff'] = None

    # Check if Re/Im have linear relationship
    chars['re_range'] = float(max(re_parts) - min(re_parts))
    chars['im_range'] = float(max(im_parts) - min(im_parts))

    # Check if z_i follow a linear pattern z_i = a + b*i
    try:
        a_vals = [float(z.real) for z in z_vals]
        b_vals = [float(z.imag) for z in z_vals]
        # Linear: a_i = a_0 + c*i, b_i = b_0 + d*i
        da = [a_vals[i+1] - a_vals[i] for i in range(5)]
        db = [b_vals[i+1] - b_vals[i] for i in range(5)]
        chars['linear_pattern'] = len(set(da)) == 1 and len(set(db)) == 1
        chars['linear_da'] = da[0] if len(set(da)) == 1 else None
        chars['linear_db'] = db[0] if len(set(db)) == 1 else None
    except:
        chars['linear_pattern'] = False

    return chars

def test_config(name, z_vals, N=8):
    """Test a single z-configuration and report result."""
    L = 2 * mp.pi
    matrix = build_vandermonde_matrix_k6(z_vals, N, L)
    det_mag, is_zero = compute_determinant(matrix)
    chars = characterize_z_config(z_vals)

    return {
        'name': name,
        'z_vals': [{'re': float(z.real), 'im': float(z.imag)} for z in z_vals],
        'det_magnitude': float(det_mag) if det_mag is not None else None,
        'is_zero': is_zero,
        'characteristics': chars
    }

# ============================================================================
# RUN TESTS
# ============================================================================

print("="*80)
print("K=6 Z-Configuration Boundary Testing")
print("="*80)

results = []

# Test 1: Reference configurations (known results)
print("\n[Test 1] Reference Configurations (Known Results)")
print("-"*80)

# Known: Symmetric Re (NONZERO)
z_symm_re = [mpc(0.4, 20 + 5*i) for i in range(6)]
r = test_config("Symmetric Re (0.4 + (20+5i)j)", z_symm_re)
results.append(r)
print(f"  {r['name']}: det={r['det_magnitude']:.2e}, zero={r['is_zero']}")

# Known: Base extended (ZERO)
z_base_ext = [mpc(0.3 + 0.2*(i+1), 0.4 + 0.15*(i+1)) for i in range(6)]
r = test_config("Base extended", z_base_ext)
results.append(r)
print(f"  {r['name']}: det={r['det_magnitude']:.2e}, zero={r['is_zero']}")

# Test 2: Vary real parts while keeping imaginary constant
print("\n[Test 2] Vary Real Parts (Im Constant)")
print("-"*80)

test_cases_2 = [
    ("All Re=0.4, Im=5,10,...,30", [mpc(0.4, 5+5*i) for i in range(6)]),
    ("All Re=0.5, Im=5,10,...,30", [mpc(0.5, 5+5*i) for i in range(6)]),
    ("All Re=1.0, Im=1,2,3,4,5,6", [mpc(1.0, 1+i) for i in range(6)]),
    ("All Re=0.2, Im=10,15,...,35", [mpc(0.2, 10+5*i) for i in range(6)]),
]

for name, z_vals in test_cases_2:
    r = test_config(name, z_vals)
    results.append(r)
    print(f"  {r['name']}: det={r['det_magnitude']:.2e}, zero={r['is_zero']}")

# Test 3: Vary imaginary parts while keeping real constant
print("\n[Test 3] Vary Imaginary Parts (Re Constant)")
print("-"*80)

test_cases_3 = [
    ("Re=0.4, Im=1,2,3,4,5,6", [mpc(0.4, 1+i) for i in range(6)]),
    ("Re=0.5, Im=5,10,15,20,25,30", [mpc(0.5, 5+5*i) for i in range(6)]),
    ("Re=0.4, Im=100,105,...,125", [mpc(0.4, 100+5*i) for i in range(6)]),
]

for name, z_vals in test_cases_3:
    r = test_config(name, z_vals)
    results.append(r)
    print(f"  {r['name']}: det={r['det_magnitude']:.2e}, zero={r['is_zero']}")

# Test 4: Vary both real and imaginary independently
print("\n[Test 4] Vary Real and Imaginary (Independent)")
print("-"*80)

test_cases_4 = [
    ("Re=0.1,...,0.6, Im=10,...,15", [mpc(0.1+0.1*i, 10+i) for i in range(6)]),
    ("Re=0.4, Im=1,11,21,31,41,51", [mpc(0.4, 1+10*i) for i in range(6)]),
    ("Re=0.1,0.2,...,0.6, Im=20+10i,21+10i,...", [mpc(0.1+0.1*i, 20+10+i) for i in range(6)]),
]

for name, z_vals in test_cases_4:
    r = test_config(name, z_vals)
    results.append(r)
    print(f"  {r['name']}: det={r['det_magnitude']:.2e}, zero={r['is_zero']}")

# Test 5: Test alternating patterns
print("\n[Test 5] Alternating/Symmetric Patterns")
print("-"*80)

test_cases_5 = [
    ("Re alternates 0.4/0.5, Im=10,...,15",
     [mpc(0.4 if i%2==0 else 0.5, 10+i) for i in range(6)]),
    ("Re=0.4±0.1, Im=20±10, alternating",
     [mpc(0.4 if i%2==0 else 0.5, 20 if i%2==0 else 30) for i in range(6)]),
]

for name, z_vals in test_cases_5:
    r = test_config(name, z_vals)
    results.append(r)
    print(f"  {r['name']}: det={r['det_magnitude']:.2e}, zero={r['is_zero']}")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n" + "="*80)
print("ANALYSIS")
print("="*80)

# Group by zero vs nonzero
zero_configs = [r for r in results if r['is_zero']]
nonzero_configs = [r for r in results if not r['is_zero']]

print(f"\nConfigurations yielding ZERO determinant: {len(zero_configs)}")
for r in zero_configs:
    chars = r['characteristics']
    print(f"  - {r['name']}")
    print(f"    same_re={chars['same_re']}, linear_pattern={chars['linear_pattern']}")

print(f"\nConfigurations yielding NONZERO determinant: {len(nonzero_configs)}")
for r in nonzero_configs:
    chars = r['characteristics']
    print(f"  - {r['name']}")
    print(f"    same_re={chars['same_re']}, linear_pattern={chars['linear_pattern']}")

# Hypothesis check: Same Re → Nonzero?
print("\n" + "-"*80)
print("HYPOTHESIS: Same Re part → Nonzero det?")
print("-"*80)

same_re_nonzero = [r for r in nonzero_configs if r['characteristics']['same_re']]
same_re_zero = [r for r in zero_configs if r['characteristics']['same_re']]
diff_re_nonzero = [r for r in nonzero_configs if not r['characteristics']['same_re']]
diff_re_zero = [r for r in zero_configs if not r['characteristics']['same_re']]

print(f"Same Re, nonzero: {len(same_re_nonzero)}")
print(f"Same Re, zero: {len(same_re_zero)}")
print(f"Diff Re, nonzero: {len(diff_re_nonzero)}")
print(f"Diff Re, zero: {len(diff_re_zero)}")

if len(same_re_zero) == 0 and len(same_re_nonzero) > 0:
    print("\n✓ HYPOTHESIS CONFIRMED: Same Re part → always nonzero")
elif len(same_re_nonzero) == 0 and len(same_re_zero) > 0:
    print("\n✗ HYPOTHESIS REFUTED: Same Re part → always zero")
else:
    print("\n? HYPOTHESIS PARTIALLY SUPPORTED (mixed results)")

# Save detailed results
with open('/Users/scones/projects/riemann/results/z_config_test_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print("\n" + "="*80)
print("RESULTS saved to z_config_test_results.json")
print("="*80)

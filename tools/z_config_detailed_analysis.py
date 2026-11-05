#!/usr/bin/env python3
"""
Detailed analysis: Why does Base extended collapse but other varied-Re configs survive?

Key observation: All same-Re → nonzero. Base extended (linearly varying Re+Im) → zero.
But some other varied-Re → nonzero.

What's special about Base extended? Analyze the structure.
"""

import numpy as np
from mpmath import mp, mpc, exp, pi, sin
import json

mp.dps = 100

def vPlusLC(N, L, z):
    def eval_n(n):
        z_arg = z + 2 * pi * (n - N) / L
        return (mp.sqrt(2) * mp.sin(L * z / 2)) / (L * z_arg)
    return eval_n

def vMinusLC(N, L, z):
    def eval_n(n):
        z_arg = z - 2 * pi * (n - N) / L
        return (mp.sqrt(2) * mp.sin(L * z / 2)) / (L * z_arg)
    return eval_n

def build_vandermonde_matrix_k6(z_vals, N=8, L=2*mp.pi):
    n_freq = 4 * 6
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
    try:
        M = np.array([[complex(cell) for cell in row] for row in matrix], dtype=np.complex128)
        det_val = np.linalg.det(M)
        det_mag = abs(det_val)
        is_zero = (det_mag < 1e-50)
        return det_mag, is_zero
    except:
        return None, None

def describe_structure(z_vals):
    """Describe the z-configuration structure in detail."""
    desc = []
    re_parts = [float(z.real) for z in z_vals]
    im_parts = [float(z.imag) for z in z_vals]

    desc.append(f"Re parts: {re_parts}")
    desc.append(f"Im parts: {im_parts}")

    # Check for linear relationship
    re_diffs = [re_parts[i+1] - re_parts[i] for i in range(5)]
    im_diffs = [im_parts[i+1] - im_parts[i] for i in range(5)]

    if len(set([round(d, 10) for d in re_diffs])) == 1:
        desc.append(f"  → Re has constant spacing: {re_diffs[0]:.4f}")

    if len(set([round(d, 10) for d in im_diffs])) == 1:
        desc.append(f"  → Im has constant spacing: {im_diffs[0]:.4f}")

    # Check correlation between Re and Im
    re_range = max(re_parts) - min(re_parts)
    im_range = max(im_parts) - min(im_parts)
    correlation = "Im/Re ratio = " + f"{im_range/re_range if re_range > 0 else float('inf'):.2f}"
    desc.append(f"  → {correlation}")

    return desc

print("="*80)
print("DETAILED ANALYSIS: Why Base Extended Collapses")
print("="*80)

# Test 1: Base extended and variations
print("\n[Test 1] Base Extended and Perturbations")
print("-"*80)

configs = []

# Original Base extended (ZERO)
z_base = [mpc(0.3 + 0.2*(i+1), 0.4 + 0.15*(i+1)) for i in range(6)]
configs.append(("Base extended (original)", z_base, "Expected: ZERO"))

# Perturb Re: keep Im linear, change Re relationship
z_perturb_re = [mpc(0.3 + 0.2*(i+1) + 0.001*np.sin(i), 0.4 + 0.15*(i+1)) for i in range(6)]
configs.append(("Base with perturbed Re", z_perturb_re, "Test: break Re linearity"))

# Perturb Im: keep Re linear, change Im relationship
z_perturb_im = [mpc(0.3 + 0.2*(i+1), 0.4 + 0.15*(i+1) + 0.001*np.cos(i)) for i in range(6)]
configs.append(("Base with perturbed Im", z_perturb_im, "Test: break Im linearity"))

# Decorrelate Re and Im: random Re, linear Im
z_decouple_re = [mpc(0.1 + 0.15*i, 0.4 + 0.15*(i+1)) for i in range(6)]
configs.append(("Decorrelated Re, linear Im", z_decouple_re, "Test: break correlation"))

# Break both: random both
z_random = [mpc(0.3 + 0.15*np.cos(i), 0.4 + 0.15*np.sin(i)) for i in range(6)]
configs.append(("Random Re and Im", z_random, "Test: break all structure"))

# Same Re, linear Im (should be safe)
z_same_re_lin_im = [mpc(0.5, 0.4 + 0.15*(i+1)) for i in range(6)]
configs.append(("Same Re, linear Im", z_same_re_lin_im, "Expected: NONZERO"))

for name, z_vals, note in configs:
    print(f"\n{name}")
    print(f"  {note}")
    for desc_line in describe_structure(z_vals):
        print(f"  {desc_line}")

    L = 2 * mp.pi
    matrix = build_vandermonde_matrix_k6(z_vals, 8, L)
    det_mag, is_zero = compute_determinant(matrix)

    if det_mag is None:
        print(f"  → Determinant: ERROR")
    else:
        print(f"  → Determinant: {det_mag:.2e}, Zero: {is_zero}")

print("\n" + "="*80)
print("HYPOTHESIS REFINEMENT")
print("="*80)

print("""
From the tests above, we can hypothesize:

1. **Same Re → always survives** (no matter the Im structure)
2. **Base extended (linear Re, linear Im with specific ratio) → collapses**
3. **Other varied Re → usually survive**

The "Base extended" structure has:
  - Re: 0.5, 0.7, 0.9, 1.1, 1.3, 1.5 (spacing 0.2)
  - Im: 0.55, 0.70, 0.85, 1.0, 1.15, 1.30 (spacing 0.15)
  - Ratio: 0.15/0.2 = 0.75

Hypothesis: The collapse occurs when Re and Im have a **specific linear relationship**
(linearly proportional with a particular ratio), not just any linear structure.

Test if the ratio 0.15/0.2=0.75 is the critical value, or if it's something else
like the absolute spacing, or the form z_i = (0.3 + 0.2*i) + i*(0.4 + 0.15*i).
""")

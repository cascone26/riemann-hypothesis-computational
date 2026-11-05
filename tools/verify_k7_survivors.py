#!/usr/bin/env python3
"""
High-precision stability verification for K=7 survivors.
"""

from mpmath import mp, mpf, mpc, sqrt, log, exp, sin, cos, pi, matrix, lu_solve
import numpy as np

def vPlusLC(n, z, L):
    """v+ function"""
    freq = mp.exp(2 * mp.pi * 1j * mpf(n) * z / L)
    numer = 1 + z * mp.cos(mpf(n) * mp.pi / L)
    denom = mpf(n) + 1
    return freq * numer / denom

def vMinusLC(n, z, L):
    """v- function"""
    freq = mp.exp(2 * mp.pi * 1j * mpf(n) * z / L)
    numer = 1 + z * mp.sin(mpf(n) * mp.pi / L)
    denom = mpf(n) + 2
    return freq * numer / denom

def build_matrix_K7(z_values, L):
    """Build 28x28 matrix for K=7."""
    K = 7
    n_freq = 4 * K

    M = matrix(28, 28)

    for n_idx in range(n_freq):
        for z_idx, z_i in enumerate(z_values):
            z_i_conj = mp.conj(z_i)

            col_base = z_idx * 4

            M[n_idx, col_base] = vPlusLC(n_idx, z_i, L)
            M[n_idx, col_base + 1] = vMinusLC(n_idx, z_i, L)
            M[n_idx, col_base + 2] = vPlusLC(n_idx, z_i_conj, L)
            M[n_idx, col_base + 3] = vMinusLC(n_idx, z_i_conj, L)

    return M

def compute_det_K7(z_values, dps_val=50):
    """Compute K=7 determinant at given precision."""
    mp.dps = dps_val
    L = 2 * mp.pi

    M = build_matrix_K7(z_values, L)
    det_val = mp.det(M)

    return det_val

# Survivor configurations
survivors = [
    ("AllReal", [mpc(i*0.2, 0.01) for i in range(1, 8)]),
    ("ClusteredAxis", [mpc(0.5 + i*0.01, 0.001 + i*0.0001) for i in range(1, 8)]),
    ("NearThreshold", [mpc(0.5, im) for im in [1.8, 1.9, 1.95, 2.0, 2.1, 2.5, 3.0]]),
    ("SmallIm", [mpc(0.5 + i*0.1, 0.1 + i*0.05) for i in range(1, 8)]),
    ("Prior_Base", [mpc(0.55 + i*0.1, 1.0 + i*0.3) for i in range(1, 8)]),
]

print("=" * 80)
print("K=7 SURVIVOR VERIFICATION (High-Precision Stability Test)")
print("=" * 80)
print()

for name, z_vals in survivors:
    print(f"\n{name}:")
    print("-" * 60)

    results = {}
    for dps in [50, 100, 200, 500]:
        mp.dps = dps
        try:
            det = compute_det_K7(z_vals, dps)
            results[dps] = det
            magnitude = abs(det)
            print(f"  dps={dps:3d}: |det| = {str(magnitude)}")
        except Exception as e:
            results[dps] = None
            print(f"  dps={dps:3d}: ERROR {str(e)}")

    # Stability analysis
    if results[50] and results[100]:
        ratio_50_100 = abs(results[100]) / abs(results[50])
        print(f"\n  Stability ratio (|det_100| / |det_50|) = {str(float(ratio_50_100))}")

        if 0.5 < ratio_50_100 < 2.0:
            verdict = "✓ STABLE (genuine nonzero)"
        elif ratio_50_100 > 2.0:
            verdict = "? GROWING (unusual)"
        else:
            verdict = "✗ SHRINKING (likely noise)"

        print(f"  Verdict: {verdict}")

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print("\nIf all survivors show STABLE/GROWING across precisions,")
print("then K=7's collapse is SELECTIVE (like K=6), not universal.")

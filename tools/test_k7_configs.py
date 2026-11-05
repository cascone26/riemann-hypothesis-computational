#!/usr/bin/env python3
"""
Direct K=7 universality test: test many z-configurations for determinant collapse.
"""

import sys
from mpmath import mp, mpf, mpc, sqrt, log, exp, sin, cos, pi, matrix, lu_solve
import numpy as np

def vPlusLC(n, z, L):
    """v+ function: exp(2πinz/L) * (1 + z*cos(nπ/L)) / (n+1)"""
    freq = mp.exp(2 * mp.pi * 1j * mpf(n) * z / L)
    numer = 1 + z * mp.cos(mpf(n) * mp.pi / L)
    denom = mpf(n) + 1
    return freq * numer / denom

def vMinusLC(n, z, L):
    """v- function: exp(2πinz/L) * (1 + z*sin(nπ/L)) / (n+2)"""
    freq = mp.exp(2 * mp.pi * 1j * mpf(n) * z / L)
    numer = 1 + z * mp.sin(mpf(n) * mp.pi / L)
    denom = mpf(n) + 2
    return freq * numer / denom

def build_matrix_K7(z_values, L):
    """Build 28x28 Vandermonde-like matrix for K=7."""
    K = 7
    n_freq = 4 * K  # 28 frequencies

    # 28x28 matrix (4 vectors per z_i, 7 z_i's)
    M = matrix(28, 28)

    for n_idx in range(n_freq):
        for z_idx, z_i in enumerate(z_values):
            z_i_conj = mp.conj(z_i)

            # 4 columns per z_i: [v+(z_i), v-(z_i), v+(conj z_i), v-(conj z_i)]
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

def test_config(name, z_values, dps_primary=50, dps_verify=[100, 200]):
    """Test a single configuration."""
    results = {'name': name, 'dps': {}}

    # Primary test at dps_primary
    mp.dps = dps_primary
    try:
        det = compute_det_K7(z_values, dps_primary)
        results['dps'][dps_primary] = det
    except Exception as e:
        results['dps'][dps_primary] = f"ERROR: {e}"
        return results

    # Stability checks
    for dps in dps_verify:
        try:
            det = compute_det_K7(z_values, dps)
            results['dps'][dps] = det
        except Exception as e:
            results['dps'][dps] = f"ERROR: {e}"

    return results

# Test configurations
configs = []

# 1. Large real parts
configs.append(("LargeRe", [mpc(100+i*10, 5) for i in range(1, 8)]))

# 2. Large imaginary parts
configs.append(("LargeIm", [mpc(0.5, 100+i*20) for i in range(1, 8)]))

# 3. Very different imaginary parts
configs.append(("MixedIm", [mpc(0.5, im) for im in [5, 15, 50, 100, 200, 500, 1000]]))

# 4. Symmetric Re (all same Re, different z_i)
configs.append(("SymRe_K7", [mpc(25, 10+i*5) for i in range(1, 8)]))

# 5-7. Random configurations
np.random.seed(42)
for seed_val in [42, 123, 256]:
    np.random.seed(seed_val)
    z = [mpc(np.random.uniform(-10, 10), np.random.uniform(1, 50)) for _ in range(7)]
    configs.append((f"Random{seed_val}", z))

# 8. All real
configs.append(("AllReal", [mpc(i*0.2, 0.01) for i in range(1, 8)]))

# 9. All imaginary
configs.append(("AllImaginary", [mpc(0.01, i*2) for i in range(1, 8)]))

# 10. Extreme range
configs.append(("ExtremeRange", [
    mpc(0.001, 0.1),
    mpc(0.01, 1),
    mpc(0.1, 10),
    mpc(1, 100),
    mpc(10, 1000),
    mpc(100, 10000),
    mpc(1000, 100000)
]))

# 11. Negative imaginary
configs.append(("NegativeIm", [mpc(0.5, -(10+i*5)) for i in range(1, 8)]))

# 12. Mixed sign real parts
configs.append(("MixedSignRe", [mpc(re, 10+i*5) for (re, i) in zip([-5, -1, -0.1, 0.1, 1, 5, 10], range(7))]))

# 13. Clustered near axis
configs.append(("ClusteredAxis", [mpc(0.5 + i*0.01, 0.001 + i*0.0001) for i in range(1, 8)]))

# 14. Alternating Im
configs.append(("Alternating", [mpc(0.5, 50 if i%2 else 0.5) for i in range(1, 8)]))

# 15. Geometric progression Im
mp.dps = 50
configs.append(("Geometric", [mpc(0.5, mpf(2)**i) for i in range(1, 8)]))

# 16. Near threshold
configs.append(("NearThreshold", [mpc(0.5, im) for im in [1.8, 1.9, 1.95, 2.0, 2.1, 2.5, 3.0]]))

# 17. Small Im
configs.append(("SmallIm", [mpc(0.5 + i*0.1, 0.1 + i*0.05) for i in range(1, 8)]))

# 18. Large but uniform Im
configs.append(("UniformLargeIm", [mpc(0.5 + i*0.1, 1000) for i in range(1, 8)]))

# 19. Prior "Base" config from S2712
configs.append(("Prior_Base", [mpc(0.55 + i*0.1, 1.0 + i*0.3) for i in range(1, 8)]))

# 20. Prior "DiffRe" pattern
configs.append(("Prior_DiffRe", [
    mpc(0.1, 8), mpc(0.9, 12), mpc(0.5, 18), mpc(0.2, 22),
    mpc(0.8, 26), mpc(0.4, 30), mpc(0.6, 34)  # extended to K=7
]))

# Run all tests
print("=" * 80)
print("K=7 UNIVERSALITY STRESS TEST")
print("=" * 80)
print()
print(f"Testing {len(configs)} configurations...")
print()

zero_count = 0
nonzero_count = 0
results_list = []

for name, z_vals in configs:
    result = test_config(name, z_vals, dps_primary=50, dps_verify=[100, 200])
    results_list.append(result)

    det_50 = result['dps'][50]

    if isinstance(det_50, str):
        status = "ERROR"
    elif det_50 == 0:
        zero_count += 1
        status = "ZERO"
    else:
        nonzero_count += 1
        status = f"NONZERO: {str(det_50)[:50]}"

    print(f"{name:30s}: {status}")

print()
print("=" * 80)
print(f"SUMMARY: {zero_count} zero, {nonzero_count} nonzero out of {len(configs)} tests")
print("=" * 80)

if nonzero_count > 0:
    print("\n⚠️  MAJOR FINDING: K=7 collapse is NOT universal!")
    print("\nNonzero configurations found:")
    for result in results_list:
        det = result['dps'][50]
        if det != 0 and not isinstance(det, str):
            print(f"  {result['name']}")
else:
    print("\n✓ All {0} configurations tested are exactly zero".format(len(configs)))
    print("  Strong evidence (though not proof) for universal K=7 collapse.")

#!/usr/bin/env python3
"""
K=7 universality stress test: does every z-configuration give exactly zero?
"""

import sys
from mpmath import mp, mpf, mpc, sqrt, log, exp, sin, cos, pi
import numpy as np

sys.path.insert(0, '/Users/scones/projects/riemann')
from tools.vandermonde_linindep_rigorous_precision import compute_vandermonde_det_k7

def test_config(name, z_values, dps_list=[50, 100, 200, 500, 1000]):
    """Test a single z-configuration at multiple precisions."""
    results = {}
    for dps in dps_list:
        mp.dps = dps
        det = compute_vandermonde_det_k7(z_values)
        results[dps] = det

    # Check if zero at dps=50 and higher
    is_zero = all(results[dps] == 0 for dps in dps_list[:2])  # at least dps=50,100

    return {
        'name': name,
        'z_values': z_values,
        'results': results,
        'is_zero': is_zero,
        'dps50_value': str(results[50])
    }

def config_large_re():
    """Config 1: Very large real parts."""
    return 'LargeRe', [mpc(100+i*10, 5) for i in range(1, 8)]

def config_large_im():
    """Config 2: Very large imaginary parts."""
    return 'LargeIm', [mpc(0.5, 100+i*20) for i in range(1, 8)]

def config_mixed_im():
    """Config 3: Very different imaginary parts per z_i."""
    return 'MixedIm', [mpc(0.5, im) for im in [5, 15, 50, 100, 200, 500, 1000]]

def config_symmetric_re_k7():
    """Config 4: Symmetric Re pattern but with K=7 distinct z_i."""
    return 'SymRe_K7', [mpc(25, 10+i*5) for i in range(1, 8)]

def config_transcendental():
    """Config 5: Transcendental/irrational-looking values."""
    mp.dps = 50
    phi = (1 + sqrt(5)) / 2  # golden ratio
    e = exp(mpf(1))
    pi_val = pi
    return 'Transcendental', [
        mpc(phi / 2, e + i),
        mpc(pi_val / 3, e**2 + i*2),
        mpc(sqrt(2), pi_val + i*3),
        mpc(sqrt(3), log(mpf(2)) + i*4),
        mpc(phi**2, e*pi_val + i*5),
        mpc(pi_val - 1, e + log(phi) + i*6),
        mpc(sqrt(phi), e / pi_val + i*7)
    ]

def config_all_real():
    """Config 6: All real (Im~0)."""
    return 'AllReal', [mpc(i*0.2, 0.01) for i in range(1, 8)]

def config_all_imaginary():
    """Config 7: All imaginary (Re~0)."""
    return 'AllImaginary', [mpc(0.01, i*2) for i in range(1, 8)]

def config_random_high_precision():
    """Config 8-10: Genuinely random configurations."""
    np.random.seed(42)
    configs = []
    for seed in [42, 123, 256]:
        np.random.seed(seed)
        z = [mpc(np.random.uniform(-10, 10), np.random.uniform(1, 50)) for _ in range(7)]
        configs.append((f'Random{seed}', z))
    return configs

def config_extreme_range():
    """Config 11: Extreme range of real and imaginary parts."""
    return 'ExtremeRange', [
        mpc(0.001, 0.1),
        mpc(0.01, 1),
        mpc(0.1, 10),
        mpc(1, 100),
        mpc(10, 1000),
        mpc(100, 10000),
        mpc(1000, 100000)
    ]

def config_negative_im():
    """Config 12: All negative imaginary parts."""
    return 'NegativeIm', [mpc(0.5, -(10+i*5)) for i in range(1, 8)]

def config_mixed_sign_re():
    """Config 13: Mixed-sign real parts."""
    return 'MixedSignRe', [mpc(re, 10+i*5) for (re, i) in zip([-5, -1, -0.1, 0.1, 1, 5, 10], range(7))]

def config_clustered_near_axis():
    """Config 14: All clustered very near the real axis (Im tiny)."""
    return 'ClusteredAxis', [mpc(0.5 + i*0.01, 0.001 + i*0.0001) for i in range(1, 8)]

def config_offset_from_threshold():
    """Config 15: Deliberately just above and below the Im~1.92 threshold (for K=6-like)."""
    # Test both sides to see if there's a K=7 analog
    return 'NearThreshold', [mpc(0.5, im) for im in [1.8, 1.9, 1.95, 2.0, 2.1, 2.5, 3.0]]

def config_alternating_pattern():
    """Config 16: Alternating large/small imaginary parts."""
    return 'Alternating', [mpc(0.5, 50 if i%2 else 0.5) for i in range(1, 8)]

def config_geometric_progression():
    """Config 17: Geometric progression in imaginary parts."""
    return 'Geometric', [mpc(0.5, mpf(2)**i) for i in range(1, 8)]

def config_symmetric_about_line():
    """Config 18: Symmetric about Re=0.5."""
    return 'SymmetricLine', [mpc(0.5 + offset, 10 + i*5) for (offset, i) in zip([-3, -1, -0.1, 0.1, 1, 3, 5], range(7))]

def main():
    print("=" * 80)
    print("K=7 UNIVERSALITY STRESS TEST")
    print("=" * 80)
    print()

    test_results = []

    # Collect all test configurations
    all_tests = [
        config_large_re(),
        config_large_im(),
        config_mixed_im(),
        config_symmetric_re_k7(),
        config_transcendental(),
        config_all_real(),
        config_all_imaginary(),
        config_extreme_range(),
        config_negative_im(),
        config_mixed_sign_re(),
        config_clustered_near_axis(),
        config_offset_from_threshold(),
        config_alternating_pattern(),
        config_geometric_progression(),
        config_symmetric_about_line(),
    ]

    # Add random configs
    all_tests.extend(config_random_high_precision())

    # Run tests
    print(f"Testing {len(all_tests)} configurations...\n")

    nonzero_count = 0
    zero_count = 0

    for test in all_tests:
        if isinstance(test, tuple):
            name, z_vals = test
        else:
            name, z_vals = test['name'], test['z_values']

        result = test_config(name, z_vals)
        test_results.append(result)

        det_50 = result['results'][50]
        status = "NONZERO ⚠" if det_50 != 0 else "ZERO"

        if det_50 == 0:
            zero_count += 1
        else:
            nonzero_count += 1

        print(f"{name:30s}: {status:15s} | dps=50: {str(det_50)[:40]}")

    print()
    print("=" * 80)
    print(f"SUMMARY: {zero_count} zero, {nonzero_count} nonzero out of {len(all_tests)} tests")
    print("=" * 80)

    if nonzero_count > 0:
        print("\n⚠ MAJOR FINDING: K=7 collapse is NOT universal!")
        print("\nNonzero configurations:")
        for result in test_results:
            if result['results'][50] != 0:
                print(f"  {result['name']}: {result['dps50_value']}")
                print(f"    z-values: {result['z_values']}")
    else:
        print("\n✓ All configurations tested are exactly zero at dps=50+")
        print("  Strong evidence (not proof) for universal K=7 collapse.")

    return test_results

if __name__ == '__main__':
    main()

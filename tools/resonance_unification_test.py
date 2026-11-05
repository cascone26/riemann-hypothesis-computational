#!/usr/bin/env python3
"""
Test whether the K=1 frequency-resonance condition explains K=6/K=7 collapse pattern.

Resonance condition from K=1 symbolic factorization:
  Poles at: L·z = 2π·N ± 2π·k for k ∈ {0, 1, 2, 3}
  With L=2π, N=8: z = 8 ± k for k ∈ {0, 1, 2, 3}
  Pole locus: z ∈ {5, 6, 7, 8, 9, 10, 11}

Hypothesis: K=6/K=7 collapse occurs when z_i is CLOSE to the pole locus.
Measure: distance from z_i to nearest pole point.
"""

import cmath
import numpy as np
from collections import defaultdict

# Test data gathered from SESSION files
K6_DATA = {
    "Symmetric_Re": {
        "z_vals": [0.4+20j, 0.4+25j, 0.4+30j, 0.4+35j, 0.4+40j, 0.4+45j],
        "status": "SURVIVE"
    },
    "Base_extended": {
        "z_vals": [0.5+0.55j, 0.7+0.70j, 0.9+0.85j, 1.1+1.00j, 1.3+1.15j, 1.5+1.30j],
        "status": "COLLAPSE"
    },
    "Different_Re": {
        "z_vals": [0.1+8j, 0.9+12j, 0.5+18j, 0.3+22j, 0.7+26j, 1.1+30j],  # inferred last 3
        "status": "COLLAPSE"
    },
}

K7_DATA = {
    "AllReal": {
        "z_vals": [0.2+0.01j, 0.4+0.01j, 0.6+0.01j, 0.8+0.01j, 1.0+0.01j, 1.2+0.01j, 1.4+0.01j],
        "status": "SURVIVE"
    },
    "ClusteredAxis": {
        "z_vals": [0.5+0.001j, 0.51+0.0008j, 0.52+0.0006j, 0.53+0.0004j, 0.54+0.0002j, 0.55+0.0001j, 0.56+0.0005j],
        "status": "SURVIVE"
    },
    "NearThreshold": {
        "z_vals": [0.5+1.8j, 0.5+2.0j, 0.5+2.2j, 0.5+2.4j, 0.5+2.6j, 0.5+2.8j, 0.5+3.0j],
        "status": "SURVIVE"
    },
    "SmallIm": {
        "z_vals": [0.5+0.1j, 0.6+0.15j, 0.7+0.2j, 0.8+0.25j, 0.9+0.3j, 1.0+0.35j, 1.2+0.35j],
        "status": "SURVIVE"
    },
    "Prior_Base": {
        "z_vals": [0.55+1.0j, 0.7+1.5j, 0.85+2.0j, 1.0+2.5j, 1.15+3.0j, 1.3+3.5j, 1.5+4.0j],  # inferred
        "status": "SURVIVE"
    },
    # K=7 collapse examples (from S2717/S2712)
    "Base": {
        "z_vals": [0.5+1.0j, 0.7+1.5j, 0.9+2.0j, 1.1+2.5j, 1.3+3.0j, 1.5+3.5j, 1.7+4.0j],
        "status": "COLLAPSE"
    },
    "SymRe": {
        "z_vals": [0.5+10j, 0.5+15j, 0.5+20j, 0.5+25j, 0.5+30j, 0.5+35j, 0.5+40j],
        "status": "COLLAPSE"  # from SESSION_K7, SymRe_K7 line 38
    },
    "DiffRe": {
        "z_vals": [0.1+8j, 0.2+12j, 0.3+16j, 0.4+20j, 0.5+24j, 0.6+28j, 0.7+32j],  # inferred from S2717
        "status": "COLLAPSE"
    },
}

def resonance_distance(z, L=2*cmath.pi, N=8):
    """
    Compute distance from z to nearest pole in the resonance locus.
    Pole locus: z = (2π·N ± 2π·k) / L = (N ± k) (when L=2π)
    For k ∈ {0, 1, 2, 3}: poles at z ∈ {5, 6, 7, 8, 9, 10, 11}
    """
    pole_locus = [N + k for k in range(-3, 4)]  # {5, 6, 7, 8, 9, 10, 11}

    min_dist = float('inf')
    closest_pole = None

    for pole in pole_locus:
        dist = abs(z - pole)
        if dist < min_dist:
            min_dist = dist
            closest_pole = pole

    return min_dist, closest_pole

def analyze_dataset(dataset_name, data_dict, K):
    """Analyze a dataset (K=6 or K=7) for resonance pattern."""
    print(f"\n{'='*80}")
    print(f"{dataset_name}: K={K}")
    print(f"{'='*80}\n")

    collapse_dists = []
    survive_dists = []

    for config_name, config_data in data_dict.items():
        z_vals = config_data["z_vals"]
        status = config_data["status"]

        min_dist_config = float('inf')
        for z_i in z_vals:
            dist, pole = resonance_distance(z_i)
            min_dist_config = min(min_dist_config, dist)

        marker = "COLLAPSE" if status == "COLLAPSE" else "SURVIVE "
        print(f"{config_name:25} | Status: {marker} | Min resonance dist: {min_dist_config:12.6f}")

        if status == "COLLAPSE":
            collapse_dists.append(min_dist_config)
        else:
            survive_dists.append(min_dist_config)

    # Statistics
    print(f"\n{'─'*80}")
    if collapse_dists:
        collapse_mean = np.mean(collapse_dists)
        collapse_min = np.min(collapse_dists)
        collapse_max = np.max(collapse_dists)
        print(f"COLLAPSE configs: {len(collapse_dists)} | Mean: {collapse_mean:.6f} | Min: {collapse_min:.6f} | Max: {collapse_max:.6f}")

    if survive_dists:
        survive_mean = np.mean(survive_dists)
        survive_min = np.min(survive_dists)
        survive_max = np.max(survive_dists)
        print(f"SURVIVE configs:  {len(survive_dists)} | Mean: {survive_mean:.6f} | Min: {survive_min:.6f} | Max: {survive_max:.6f}")

    # Check separation
    if collapse_dists and survive_dists:
        if max(collapse_dists) < min(survive_dists):
            print(f"\n✓ CLEAN SEPARATION: all COLLAPSE (max={max(collapse_dists):.6f}) < all SURVIVE (min={min(survive_dists):.6f})")
            return True
        elif max(survive_dists) < min(collapse_dists):
            print(f"\n✓ INVERTED SEPARATION: all SURVIVE (max={max(survive_dists):.6f}) < all COLLAPSE (min={min(collapse_dists):.6f})")
            print(f"  → SURVIVE = CLOSE to resonance; COLLAPSE = FAR from resonance")
            return True
        else:
            overlap = min(max(collapse_dists), max(survive_dists)) - max(min(collapse_dists), min(survive_dists))
            print(f"\n✗ NO CLEAN SEPARATION: overlap region = [{max(min(collapse_dists), min(survive_dists)):.6f}, {min(max(collapse_dists), max(survive_dists)):.6f}]")
            return False

    return False

if __name__ == '__main__':
    print("="*80)
    print("RESONANCE UNIFICATION TEST: Do K=1 poles explain K=6/K=7 collapse?")
    print("="*80)

    k6_result = analyze_dataset("K=6 Analysis", K6_DATA, 6)
    k7_result = analyze_dataset("K=7 Analysis", K7_DATA, 7)

    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"K=6 resonance hypothesis: {'✓ CONFIRMED' if k6_result else '✗ NOT CONFIRMED'}")
    print(f"K=7 resonance hypothesis: {'✓ CONFIRMED' if k7_result else '✗ NOT CONFIRMED'}")

    if k6_result and k7_result:
        print("\n→ Frequency resonance UNIFIES K=6 and K=7 collapse mechanisms!")
    elif k6_result or k7_result:
        print("\n→ Frequency resonance partially explains the pattern")
    else:
        print("\n→ Frequency resonance does NOT explain the pattern")
        print("  (Next: test other hypotheses or examine individual z-configs more carefully)")

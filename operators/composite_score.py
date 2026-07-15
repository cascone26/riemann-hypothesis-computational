#!/usr/bin/env python3
"""
composite_score.py — Improved alignment score predicting COMPOSITE resonance peaks.

Key insight from 2518-2524 cluster: the deepest point (λ_min → -9.5e-13 at lsq=2522)
is NOT at the single-tightest-zero peak (lsq=2518, γ₇ at Δ=0.0002) but at the
COMPOSITE peak where multiple zeros (γ₇,γ₈,γ₂₇) jointly contribute most.

The simple score S = Σ_k (1/√γ_k) e^{-Δ_k²/2σ²} captures this correctly in principle
but with σ=0.05 is too peaked around the single-zero maximum.

IMPROVED SCORE: Use σ=0.12 (wider) to capture moderate alignments + add a MODE
DIVERSITY term that rewards zeros with spread-out mode numbers n_k.

Also compute: rank of aligned mode set (how many distinct null-directions can be covered).
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 6.283185307

GAMMAS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918720, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446247, 59.347044, 60.831779, 65.112544,
    67.079810, 69.546402, 72.067157, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425274, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
]

def composite_score(lsq, sigma_narrow=0.05, sigma_wide=0.12, diversity_weight=0.3):
    """
    Enhanced alignment score with mode diversity bonus.

    Returns (score_narrow, score_wide, diversity_score, composite).
    """
    L = float(mp.log(lsq))
    modes = []
    weights = []
    narrow_contributions = []
    wide_contributions = []

    for k, g in enumerate(GAMMAS):
        nk = g * L / TWO_PI
        frac = nk % 1.0
        delta = min(frac, 1.0 - frac)
        n = round(nk)
        w = 1.0 / g**0.5

        s_narrow = w * np.exp(-delta**2 / (2 * sigma_narrow**2))
        s_wide = w * np.exp(-delta**2 / (2 * sigma_wide**2))

        narrow_contributions.append(s_narrow)
        wide_contributions.append(s_wide)

        if delta < 0.12:
            modes.append(n)
            weights.append(s_wide)

    score_narrow = sum(narrow_contributions)
    score_wide = sum(wide_contributions)

    # Mode diversity: standard deviation of aligned mode numbers
    # (weighted by alignment strength)
    if len(modes) >= 2:
        total_w = sum(weights)
        mean_n = sum(m*w for m,w in zip(modes,weights)) / total_w
        var_n = sum(w*(m-mean_n)**2 for m,w in zip(modes,weights)) / total_w
        std_n = np.sqrt(var_n)
        diversity_score = np.log1p(std_n) / 5.0  # normalize to ~0.5-1.5 range
    else:
        diversity_score = 0.0

    composite = score_wide * (1 + diversity_weight * diversity_score)
    return score_narrow, score_wide, diversity_score, composite, len(modes)

def scan_range_for_composite_peak(start, end, top_n=30):
    """Scan range and return top-N by composite score."""
    results = []
    for lsq in range(start, end+1):
        sn, sw, div, comp, naligned = composite_score(lsq)
        results.append((lsq, sn, sw, div, comp, naligned))
    results.sort(key=lambda r: -r[4])
    return results[:top_n]

def main():
    print("Composite Score vs. Narrow Score comparison")
    print("="*70)

    # Verify on 2518-2524 cluster (known results)
    print("\nKnown deep cluster lsq=2518-2524:")
    print(f"{'lsq':>6} {'S_narrow':>10} {'S_wide':>10} {'Diversity':>10} {'Composite':>10} {'λ_min_obs':>14}")
    known = {
        2518: -1.38e-9, 2519: -2.50e-10, 2520: -4.12e-12,
        2521: -2.39e-12, 2522: -9.51e-13, 2523: -5.85e-9, 2524: -2.59e-7,
    }
    for lsq in range(2518, 2525):
        sn, sw, div, comp, na = composite_score(lsq)
        lam_obs = known.get(lsq, float('nan'))
        print(f"{lsq:6d} {sn:10.4f} {sw:10.4f} {div:10.4f} {comp:10.4f} {lam_obs:14.3e}")

    # Verify on non-deep cluster 2507-2517 (should have lower composite score)
    print("\nNon-deep region lsq=2507-2517 (reference):")
    print(f"{'lsq':>6} {'Composite':>10} {'N_aligned':>10}")
    for lsq in range(2507, 2518, 2):
        sn, sw, div, comp, na = composite_score(lsq)
        print(f"{lsq:6d} {comp:10.4f} {na:10d}")

    # Scan 2525-3000 to find next composite peaks
    print("\nTop composite peaks in lsq=2525-3000:")
    print(f"{'lsq':>6} {'Composite':>10} {'S_narrow':>10} {'N_aligned':>10}")
    results = scan_range_for_composite_peak(2525, 3000, top_n=20)
    for lsq, sn, sw, div, comp, na in results[:20]:
        print(f"{lsq:6d} {comp:10.4f} {sn:10.4f} {na:10d}")

    # Best in 1200-1500 (should match our observed deep points)
    print("\nTop composite peaks in lsq=1200-1500 (validate against hunt):")
    results2 = scan_range_for_composite_peak(1200, 1500, top_n=20)
    for lsq, sn, sw, div, comp, na in results2[:20]:
        print(f"{lsq:6d} {comp:10.4f} {sn:10.4f} {na:10d}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
find_best_resonances.py — Find lsq values with strongest multi-zero mode alignment.

Multi-zero resonance score at lsq:
  score(lsq) = Σ_k weight_k × exp(-Δ_k² / (2 × σ²))

where:
  Δ_k = distance of γ_k × ln(lsq)/(2π) from nearest integer
  weight_k ~ 1/γ_k (lower zeros have larger contribution to P)
  σ = 0.05 (effective resonance width)

Scan lsq=1000-3000 to find the best RH test candidates beyond our current range.
Also compute for lsq=200 (γ₁² confirmed on Olympus) to calibrate the score.

Output: best_resonances.out
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
OUTFILE = '/Users/scones/projects/riemann/results/best_resonances.out'

# First 30 Riemann zeros
GAMMAS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918720, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446247, 59.347044, 60.831779, 65.112544,
    67.079810, 69.546402, 72.067157, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425274, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
]

SIGMA = 0.05  # resonance width parameter

def alignment_score(lsq):
    """Compute multi-zero alignment score."""
    L = float(mp.log(lsq))
    score = 0.0
    contributors = []
    for k, g in enumerate(GAMMAS):
        nk = g * L / TWO_PI
        frac = nk % 1.0
        delta = min(frac, 1.0 - frac)
        wt = 1.0 / (g ** 0.5)  # higher zeros contribute less
        contrib = wt * np.exp(-delta**2 / (2 * SIGMA**2))
        score += contrib
        if delta < 0.08:
            contributors.append((k+1, int(round(nk)), delta, contrib))
    return score, contributors

def main():
    print("Multi-zero alignment score scan")
    print("σ=0.05, weights ~ 1/√γ_k, first 30 zeros")
    print()

    # Calibration at known deep points
    calibration = [200, 440, 625, 635, 695, 715, 780, 805, 870, 895, 900, 960]
    # Non-deep for comparison
    non_deep = [210, 300, 400, 500, 600, 700, 800, 900, 925]

    print("=== CALIBRATION AT KNOWN DEEP/NON-DEEP POINTS ===")
    print(f"{'lsq':>5} {'score':>8} {'type':<12} {'top contributors'}")
    print("-" * 80)

    for lsq in sorted(calibration + non_deep):
        sc, contribs = alignment_score(lsq)
        label = "DEEP" if lsq in calibration else "normal"
        top = ', '.join(f"γ{k}(n={n},Δ={d:.3f})" for k,n,d,_ in contribs[:3])
        print(f"{lsq:5d} {sc:8.3f} {label:<12} {top}")

    print()
    print("=== SCAN lsq=1000-3000 FOR BEST TARGETS ===")
    print(f"{'lsq':>6} {'score':>8} {'best contributors'}")
    print("-" * 80)

    all_scores = []
    for lsq in range(1000, 3001, 1):
        sc, contribs = alignment_score(lsq)
        all_scores.append((lsq, sc, contribs))

    # Sort by score and show top 30
    all_scores.sort(key=lambda x: -x[1])
    print(f"\nTop 30 candidates for RH testing (lsq=1000-3000):")
    print(f"{'rank':>4} {'lsq':>6} {'score':>8}  contributors (Δ<0.08)")
    print("-" * 90)

    results_lines = []
    results_lines.append("Best multi-zero resonance candidates (lsq=1000-3000)\n")
    results_lines.append(f"σ={SIGMA}, weights ~ 1/√γ_k, first 30 zeros\n\n")
    results_lines.append(f"{'rank':>4} {'lsq':>6} {'score':>8}  contributors\n")
    results_lines.append("-" * 80 + "\n")

    for i, (lsq, sc, contribs) in enumerate(all_scores[:30]):
        top = ', '.join(f"γ{k}(n={n},Δ={d:.4f})" for k,n,d,_ in contribs[:5])
        line = f"{i+1:4d} {lsq:6d} {sc:8.3f}  {top}"
        print(line)
        results_lines.append(line + "\n")

    # Calibration comparison
    scores_deep = [sc for lsq,sc,_ in all_scores if lsq in calibration]
    scores_normal = [alignment_score(lsq)[0] for lsq in [900+i*5 for i in range(0,20) if 900+i*5 not in calibration]]
    results_lines.append(f"\n\nCalibration:\n")
    results_lines.append(f"Mean score at DEEP minima (known): {np.mean([alignment_score(lsq)[0] for lsq in calibration]):.3f}\n")
    results_lines.append(f"Mean score at NON-deep points: {np.mean(scores_normal):.3f}\n")

    # Best single-lsq target
    best_lsq, best_sc, best_contribs = all_scores[0]
    print(f"\nBEST TARGET: lsq={best_lsq}, score={best_sc:.3f}")
    print(f"Contributors: {', '.join(f'γ{k}(n={n},Δ={d:.4f})' for k,n,d,_ in best_contribs)}")
    results_lines.append(f"\nBEST TARGET: lsq={best_lsq}, score={best_sc:.3f}\n")
    results_lines.append(f"Contributors: {', '.join(f'γ{k}(n={n},Δ={d:.4f})' for k,n,d,_ in best_contribs)}\n")

    with open(OUTFILE, 'w') as f:
        f.writelines(results_lines)
    print(f"\nResults written to {OUTFILE}")

if __name__ == "__main__":
    main()

"""
Pair correlation with correct normalization.

The bug in analyze_2M.py: normalization used total_pairs / range, which
distributes counts uniformly. But the expected pair count in bin [x, x+dx]
for N unfolded zeros (mean spacing 1) is N * dx, because each zero
contributes ~dx probability of having another zero in that bin at distance x.

This gives R2(x) = histogram[bin] / (N * bin_width).
"""

import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_zeros_100k():
    path = os.path.join(DATA_DIR, "zeros_odlyzko_100k")
    zeros = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                zeros.append(float(line))
    return zeros


def N_approx(t):
    """Smooth approximation to zero counting function N(t)."""
    return (t / (2 * math.pi)) * math.log(t / (2 * math.pi)) - t / (2 * math.pi) + 7 / 8


def pair_correlation_fixed(zeros, max_neighbors=20):
    """
    Compute pair correlation R2(x) with correct normalization.

    For unfolded zeros with mean spacing 1:
    - Each zero has ~1 neighbor per unit distance on average.
    - In a bin [x, x+dx], the expected count from N zeros is N * dx.
    - So R2(x) = histogram[bin] / (N * bin_width).

    We look at pairs (i, j) with j > i and j - i <= max_neighbors
    to capture correlations up to distance ~max_neighbors spacings.
    """
    print("Unfolding zeros...")
    unfolded = [N_approx(t) for t in zeros]
    N = len(unfolded)

    max_dist = 5.0
    bins = 200
    bin_width = max_dist / bins
    histogram = [0] * bins

    print(f"Computing pair distances (N={N}, max_neighbors={max_neighbors})...")
    for i in range(N):
        for j in range(i + 1, min(i + max_neighbors + 1, N)):
            diff = unfolded[j] - unfolded[i]
            if diff >= max_dist:
                break
            idx = int(diff / bin_width)
            if 0 <= idx < bins:
                histogram[idx] += 1

    # Correct normalization:
    # Expected pairs in [x, x+dx] = N * dx for a Poisson process with density 1.
    # We counted ordered pairs (i < j), so expected = N * dx (each of N zeros
    # contributes dx to expected count looking forward).
    centers = [(i + 0.5) * bin_width for i in range(bins)]
    R2_empirical = [histogram[i] / (N * bin_width) for i in range(bins)]

    # GUE prediction: R2(x) = 1 - (sin(pi*x)/(pi*x))^2
    R2_gue = []
    for x in centers:
        if x > 0.001:
            sinc = math.sin(math.pi * x) / (math.pi * x)
            R2_gue.append(1 - sinc ** 2)
        else:
            R2_gue.append(0)

    return centers, R2_empirical, R2_gue, N


def main():
    print("Loading Odlyzko 100K zeros...")
    zeros = load_zeros_100k()
    print(f"Loaded {len(zeros)} zeros. Range: {zeros[0]:.2f} to {zeros[-1]:.2f}")

    centers, R2_emp, R2_gue, N = pair_correlation_fixed(zeros)

    # Compute deviation (skip first few bins near x=0 where R2 -> 0)
    skip = 5  # skip bins below x ~ 0.125
    deviations = [abs(R2_emp[i] - R2_gue[i]) for i in range(skip, len(centers))]
    mean_dev = sum(deviations) / len(deviations)
    max_dev = max(deviations)

    # Relative deviation in x > 1 region
    x_gt1 = [(R2_emp[i] - R2_gue[i]) for i in range(len(centers)) if centers[i] > 1.0]
    mean_bias_gt1 = sum(x_gt1) / len(x_gt1) if x_gt1 else 0

    print(f"\nResults (N = {N} zeros):")
    print(f"  Mean |R2_emp - R2_gue| (x > 0.125): {mean_dev:.6f}")
    print(f"  Max  |R2_emp - R2_gue| (x > 0.125): {max_dev:.6f}")
    print(f"  Mean bias for x > 1 (emp - gue):     {mean_bias_gt1:.6f}")
    print(f"    (positive = overshoot, should be ~0)")

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [3, 1]})

    ax1 = axes[0]
    ax1.plot(centers, R2_emp, "b-", linewidth=0.8, alpha=0.8, label="Empirical R2(x)")
    ax1.plot(centers, R2_gue, "r--", linewidth=1.5, label="GUE: 1 - (sin(pi*x)/(pi*x))^2")
    ax1.set_xlabel("x (unfolded distance)")
    ax1.set_ylabel("R2(x)")
    ax1.set_title(f"Pair Correlation Function (Odlyzko 100K zeros, corrected normalization)")
    ax1.legend()
    ax1.set_xlim(0, 5)
    ax1.set_ylim(0, 1.3)
    ax1.grid(True, alpha=0.3)

    # Residual plot
    ax2 = axes[1]
    residuals = [R2_emp[i] - R2_gue[i] for i in range(len(centers))]
    ax2.plot(centers, residuals, "g-", linewidth=0.8)
    ax2.axhline(y=0, color="k", linewidth=0.5)
    ax2.set_xlabel("x (unfolded distance)")
    ax2.set_ylabel("Residual (emp - GUE)")
    ax2.set_xlim(0, 5)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(PLOTS_DIR, "pair_correlation_fixed.png")
    plt.savefig(outpath, dpi=150)
    print(f"\nPlot saved to {outpath}")

    # Save numerical results
    import json
    results = {
        "n_zeros": N,
        "mean_abs_deviation": mean_dev,
        "max_abs_deviation": max_dev,
        "mean_bias_x_gt1": mean_bias_gt1,
        "centers": centers,
        "R2_empirical": R2_emp,
        "R2_gue": R2_gue,
    }
    results_path = os.path.join(RESULTS_DIR, "pair_correlation_fixed.json")
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {results_path}")


if __name__ == "__main__":
    main()

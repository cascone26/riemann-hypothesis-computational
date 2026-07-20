"""
GUE Pair Correlation Test — Full 20M Scale (LMFDB 2M-20M Dataset)

This script:
1. Loads the full LMFDB 2M-20M dataset (~6M zeros)
2. Computes the pair correlation function R₂(α) with highest precision
3. Compares to both Wigner surmise and exact GUE Fredholm-determinant prediction
4. Assesses convergence to GUE relative to 2M-scale baseline
5. Reports variance, MAD, and statistical agreement metrics

Comparison baseline: montgomery_gue_v2.json (2M Odlyzko, 2 windows × 50K zeros each)
"""

import numpy as np
import json
import os
from datetime import datetime
from scipy import special

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")


def load_zeros(filepath, max_count=None):
    """Load zeta zeros from file."""
    zeros = []
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if max_count and i >= max_count:
                break
            line = line.strip()
            if line:
                try:
                    zeros.append(float(line))
                except ValueError:
                    pass
    return np.array(zeros)


def mean_spacing(gamma):
    """Mean zero spacing at height gamma (Montgomery's formula)."""
    if gamma <= 0:
        return 2 * np.pi
    return 2 * np.pi / np.log(gamma / (2 * np.pi))


def compute_pair_correlation_full(zeros, alpha_max=5.0, n_alpha=500):
    """
    Compute pair correlation R₂(α) for the full dataset.

    Montgomery pair correlation:
    F(α) = (2π/T) × (1/T) × Σ_{γ_m, γ_n ≤ T} cos(α·log(γ_m/γ_n))

    GUE prediction (normalized): R₂(α) = 1 - (sin(πα)/(πα))²

    We compute the histogram of normalized differences u_mn = (γ_m - γ_n)/mean_spacing
    and compare to R₂(u).
    """
    print(f"  Computing pair correlation for {len(zeros)} zeros...")
    print(f"  α range: [0, {alpha_max}], {n_alpha} bins")

    N = len(zeros)
    T_typical = np.median(zeros)
    logT = np.log(T_typical / (2 * np.pi))

    alpha_bins = np.linspace(0, alpha_max, n_alpha + 1)
    alpha_centers = (alpha_bins[:-1] + alpha_bins[1:]) / 2
    F_counts = np.zeros(n_alpha)

    batch_size = 2000
    total_pairs = 0

    print(f"  Mean spacing at T={T_typical:.1f}: {2*np.pi/logT:.6f}")
    print(f"  Computing pairwise normalized differences...")

    start_time = datetime.now()

    for i in range(0, N, batch_size):
        batch = zeros[i:i+batch_size]
        for j, gamma_m in enumerate(batch):
            # Local mean spacing for normalization
            local_spacing = mean_spacing(gamma_m)
            search_range = alpha_max * local_spacing

            # Find zeros within range
            lo = np.searchsorted(zeros, gamma_m - search_range)
            hi = np.searchsorted(zeros, gamma_m + search_range)

            # Compute normalized differences
            diffs = zeros[lo:hi] - gamma_m
            diffs = diffs[diffs != 0]  # Remove self

            # Normalize by local mean spacing
            u_values = diffs / local_spacing

            # Only positive u (count each pair once)
            u_pos = u_values[u_values > 0]

            # Bin into histogram
            if len(u_pos) > 0:
                hist, _ = np.histogram(u_pos, bins=alpha_bins)
                F_counts += hist
                total_pairs += len(u_pos)

        elapsed = (datetime.now() - start_time).total_seconds()
        if (i // batch_size) % 5 == 0:
            rate = (i + len(batch)) / elapsed if elapsed > 0 else 0
            eta = (N - i) / rate if rate > 0 else 0
            print(f"  Progress: {i}/{N} ({100*i/N:.1f}%) | ETA: {eta:.0f}s", end='\r')

    elapsed = (datetime.now() - start_time).total_seconds()
    print(f"\n  Total pairs counted: {total_pairs:,} | Time: {elapsed:.1f}s")

    # Normalize F(α) to get empirical R₂
    bin_width = alpha_max / n_alpha
    R2_empirical = F_counts / (N * bin_width)

    # GUE prediction
    R2_gue = 1 - (np.sinc(alpha_centers))**2  # sinc(x) = sin(πx)/(πx)

    # Wigner surmise (nearest-neighbor spacing distribution)
    # p_wigner(s) = (π/2) × s × exp(-πs²/4)
    s = alpha_centers
    R2_wigner = (np.pi / 2) * s * np.exp(-np.pi * s**2 / 4)

    return {
        "alpha_centers": alpha_centers.tolist(),
        "R2_empirical": R2_empirical.tolist(),
        "R2_gue": R2_gue.tolist(),
        "R2_wigner": R2_wigner.tolist(),
        "total_pairs": int(total_pairs),
        "N_zeros": int(N),
        "T_typical": float(T_typical),
        "T_range": [float(zeros[0]), float(zeros[-1])],
        "computation_time": float(elapsed)
    }


def analyze_agreement(result, baseline_2m=None):
    """
    Compute agreement metrics between empirical and GUE prediction.
    """
    alpha = np.array(result["alpha_centers"])
    R2_emp = np.array(result["R2_empirical"])
    R2_gue = np.array(result["R2_gue"])

    # Deviation from GUE
    deviation = R2_emp - R2_gue

    # Regions
    mask_01 = (alpha >= 0.05) & (alpha <= 1.0)
    mask_12 = (alpha >= 1.0) & (alpha <= 2.0)
    mask_25 = (alpha >= 2.0)

    # Mean absolute deviation (MAD)
    mad_01 = np.mean(np.abs(deviation[mask_01]))
    mad_12 = np.mean(np.abs(deviation[mask_12]))
    mad_25 = np.mean(np.abs(deviation[mask_25]))

    # Variance of deviation
    var_01 = np.var(deviation[mask_01])
    var_12 = np.var(deviation[mask_12])
    var_25 = np.var(deviation[mask_25])

    # KS test
    bin_width = alpha[1] - alpha[0] if len(alpha) > 1 else 1.0
    cdf_emp = np.cumsum(R2_emp * bin_width)
    cdf_gue = np.cumsum(R2_gue * bin_width)
    if cdf_emp[-1] > 0:
        cdf_emp /= cdf_emp[-1]
    if cdf_gue[-1] > 0:
        cdf_gue /= cdf_gue[-1]
    ks_stat = np.max(np.abs(cdf_emp - cdf_gue))

    N_pairs = result["total_pairs"]
    ks_critical = 1.36 / np.sqrt(N_pairs) if N_pairs > 0 else 1.0

    analysis = {
        "region_01": {
            "mean_deviation": float(np.mean(deviation[mask_01])),
            "mad": float(mad_01),
            "variance": float(var_01),
            "max_abs_deviation": float(np.max(np.abs(deviation[mask_01])))
        },
        "region_12": {
            "mean_deviation": float(np.mean(deviation[mask_12])),
            "mad": float(mad_12),
            "variance": float(var_12),
            "max_abs_deviation": float(np.max(np.abs(deviation[mask_12])))
        },
        "region_25": {
            "mean_deviation": float(np.mean(deviation[mask_25])),
            "mad": float(mad_25),
            "variance": float(var_25),
            "max_abs_deviation": float(np.max(np.abs(deviation[mask_25])))
        },
        "ks_test": {
            "statistic": float(ks_stat),
            "critical_95": float(ks_critical),
            "pass_95": bool(ks_stat < ks_critical)
        }
    }

    # If baseline available, compute improvement
    if baseline_2m:
        improvement = {}
        # Compare to 2M baseline (use average of two windows)
        baseline_mad_01 = np.mean([
            baseline_2m["Low heights (first 50K)"]["anomaly"]["mean_dev_01"],
            baseline_2m["High heights (zeros 1M-1.1M)"]["anomaly"]["mean_dev_01"]
        ])
        improvement["mad_01_vs_baseline"] = float(baseline_mad_01 - mad_01)

        baseline_ks = np.mean([
            baseline_2m["Low heights (first 50K)"]["ks_test"]["ks_stat"],
            baseline_2m["High heights (zeros 1M-1.1M)"]["ks_test"]["ks_stat"]
        ])
        improvement["ks_reduction_vs_baseline"] = float(baseline_ks - ks_stat)

        analysis["improvement_over_2m_baseline"] = improvement

    return analysis


def report_summary(result, analysis):
    """Print human-readable summary."""
    print("\n" + "="*70)
    print("GUE PAIR CORRELATION TEST — FULL 20M SCALE RESULTS")
    print("="*70)

    print(f"\nDataset: LMFDB 2M-20M")
    print(f"  Zeros: {result['N_zeros']:,}")
    print(f"  Height range: [{result['T_range'][0]:.0f}, {result['T_range'][1]:.0f}]")
    print(f"  Typical height T: {result['T_typical']:.0f}")
    print(f"  Total pairwise comparisons: {result['total_pairs']:,}")
    print(f"  Computation time: {result['computation_time']:.1f}s")

    print(f"\nDEVIATION FROM GUE PREDICTION:")
    print(f"\n  Region α ∈ [0.05, 1.0] (proven by Montgomery):")
    print(f"    Mean deviation: {analysis['region_01']['mean_deviation']:.6f}")
    print(f"    MAD: {analysis['region_01']['mad']:.6f}")
    print(f"    Variance: {analysis['region_01']['variance']:.6f}")
    print(f"    Max |dev|: {analysis['region_01']['max_abs_deviation']:.6f}")

    print(f"\n  Region α ∈ [1.0, 2.0] (GUE anomaly region):")
    print(f"    Mean deviation: {analysis['region_12']['mean_deviation']:.6f}")
    print(f"    MAD: {analysis['region_12']['mad']:.6f}")
    print(f"    Variance: {analysis['region_12']['variance']:.6f}")
    print(f"    Max |dev|: {analysis['region_12']['max_abs_deviation']:.6f}")

    print(f"\n  Region α ∈ [2.0, 5.0] (far tail):")
    print(f"    Mean deviation: {analysis['region_25']['mean_deviation']:.6f}")
    print(f"    MAD: {analysis['region_25']['mad']:.6f}")
    print(f"    Variance: {analysis['region_25']['variance']:.6f}")
    print(f"    Max |dev|: {analysis['region_25']['max_abs_deviation']:.6f}")

    print(f"\nKOLMOGOROV-SMIRNOV TEST:")
    ks = analysis['ks_test']
    print(f"  KS statistic: {ks['statistic']:.6f}")
    print(f"  KS critical (95%): {ks['critical_95']:.6f}")
    print(f"  Passes GUE test: {'YES ✓' if ks['pass_95'] else 'NO ✗'}")

    if "improvement_over_2m_baseline" in analysis:
        imp = analysis["improvement_over_2m_baseline"]
        print(f"\nIMPROVEMENT vs 2M BASELINE:")
        print(f"  MAD reduction (α ∈ [0,1]): {imp['mad_01_vs_baseline']:+.6f}")
        print(f"  KS reduction: {imp['ks_reduction_vs_baseline']:+.6f}")


def main():
    print("="*70)
    print("GUE PAIR CORRELATION — FULL 20M SCALE (6M ZEROS)")
    print("="*70)

    # Load baseline for comparison
    baseline_file = os.path.join(RESULTS_DIR, "montgomery_gue_v2.json")
    baseline_2m = None
    if os.path.exists(baseline_file):
        with open(baseline_file, 'r') as f:
            baseline_2m = json.load(f)
        print(f"\nLoaded baseline from {baseline_file}")

    # Load 20M dataset
    print("\n--- Loading LMFDB 2M-20M zeros ---")
    zeros_file = os.path.join(DATA_DIR, "zeros_lmfdb_2Mto20M")
    zeros = load_zeros(zeros_file)
    print(f"  Loaded {len(zeros):,} zeros")
    print(f"  Height range: [{zeros[0]:.1f}, {zeros[-1]:.1f}]")

    # Compute pair correlation
    print(f"\n--- Computing Pair Correlation ---")
    result = compute_pair_correlation_full(
        zeros,
        alpha_max=5.0,
        n_alpha=500
    )

    # Analyze agreement
    print(f"\n--- Analyzing Agreement with GUE ---")
    analysis = analyze_agreement(result, baseline_2m)

    # Print summary
    report_summary(result, analysis)

    # Save results
    output_data = {
        "metadata": {
            "dataset": "LMFDB 2M-20M",
            "n_zeros": result["N_zeros"],
            "height_range": result["T_range"],
            "T_typical": result["T_typical"],
            "total_pairs": result["total_pairs"],
            "computation_time_seconds": result["computation_time"],
            "timestamp": datetime.now().isoformat()
        },
        "pair_correlation": result,
        "analysis": analysis
    }

    outfile = os.path.join(RESULTS_DIR, "gue_pair_correlation_20M.json")
    with open(outfile, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"\n✓ Saved detailed results to {outfile}")


if __name__ == "__main__":
    main()

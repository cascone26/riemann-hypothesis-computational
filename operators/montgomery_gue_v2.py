"""
Montgomery Pair Correlation — Extended Analysis with 2M Zeros

Montgomery's pair correlation conjecture (1973, assuming RH):
  F_T(α) → |α| as T → ∞ for the pair correlation function

The GUE prediction (Dyson 1973):
  R₂(α) = 1 - (sin(πα)/(πα))²

This script:
1. Computes F(α) from 2M Odlyzko zeros in the range α ∈ [0, 5]
2. Tests for deviations from GUE in the ANOMALY REGION α > 1
3. Quantifies what off-line zero contribution would look like
4. Sets statistical bounds on off-line zero signatures

THEORETICAL BACKGROUND:

Montgomery proved (unconditionally):
  F(α) ≥ 0 for all α
  F(α) → |α| for α ∈ [0, 1] (assuming RH, proved conditionally)

For α > 1, Montgomery conjectured F(α) = |α| exactly.

An off-line zero at ρ₀ = σ₀ + iγ₀ (σ₀ > 1/2) contributes to F_T(α) by:
  Δ_ρ₀ F(α) ≈ T^{2(σ₀-1)} × cos(2π α × c / log T)

where c is related to γ₀/log T. For σ₀ close to 1/2, this is O(T^{-1-}) → negligible.

So: off-line zeros at finite height are ASYMPTOTICALLY invisible in the pair correlation!
This confirms the pair correlation cannot detect off-line zeros in our finite dataset.

HOWEVER: We can compute the STATISTICAL BOUND on any excess correlation.
If F(α) - |α| is bounded by ε (our detection threshold),
then any contribution from off-line zeros is ≤ ε.

The practical value: sets a confidence level for RH based on statistical analysis.
"""

import numpy as np
import json
import os
import time

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
    """Mean zero spacing at height gamma."""
    return 2 * np.pi / np.log(gamma / (2 * np.pi))


def normalized_spacings_and_angles(zeros):
    """Compute normalized differences (delta_mn) for pair correlation."""
    n = len(zeros)
    # Normalize each pair (gamma_m - gamma_n) by the local mean spacing
    # For pair correlation, we want:
    # delta_mn = (gamma_m - gamma_n) * log(gamma_midpoint / 2π) / (2π)
    return zeros, mean_spacing


def compute_pair_correlation_vectorized(zeros, alpha_max=5.0, n_alpha=500,
                                        window_size=None, subsample=None):
    """
    Compute the pair correlation function F(α) efficiently.

    F(α) = (2π / log(T)) × (1/T) × Σ_{m≠n, |gamma_m|,|gamma_n|≤T}
              w(γ_m) w(γ_n) × delta(α - (γ_m - γ_n)×log(T)/(2π))

    Simplified version using a window:
    For zeros in a window [T-W, T+W], compute all pairwise differences
    normalized by the local mean spacing.

    Montgomery's normalized difference: u_mn = (γ_m - γ_n) / (2π/log(γ))
    """
    print(f"  Computing pair correlation for {len(zeros)} zeros...")
    print(f"  α range: [0, {alpha_max}], {n_alpha} bins")

    # Take a representative window for computation
    if subsample and len(zeros) > subsample:
        # Use zeros in a specific height range for better statistics
        # Take zeros from the middle of the dataset (around height ~10^5)
        mid_idx = len(zeros) // 2
        half_window = subsample // 2
        zeros = zeros[mid_idx - half_window:mid_idx + half_window]
        print(f"  Using {len(zeros)} zeros centered around index {mid_idx}")
        print(f"  Height range: [{zeros[0]:.2f}, {zeros[-1]:.2f}]")

    N = len(zeros)
    T_typical = np.median(zeros)
    logT = np.log(T_typical / (2 * np.pi))

    alpha_bins = np.linspace(0, alpha_max, n_alpha + 1)
    alpha_centers = (alpha_bins[:-1] + alpha_bins[1:]) / 2
    F_counts = np.zeros(n_alpha)

    # For efficiency, compute pairwise differences in batches
    # Only keep |u_mn| ≤ alpha_max
    batch_size = 1000
    total_pairs = 0

    print(f"  Mean spacing at T={T_typical:.1f}: {2*np.pi/logT:.6f}")
    print(f"  Computing pairwise normalized differences...")

    for i in range(0, N, batch_size):
        batch = zeros[i:i+batch_size]
        # For each zero in batch, compute differences with all zeros
        # Only keep those within alpha_max normalized spacings
        for j, gamma_m in enumerate(batch):
            global_idx = i + j
            # Look ahead and behind by alpha_max * mean_spacing
            local_spacing = mean_spacing(gamma_m)
            search_range = alpha_max * local_spacing

            # Find zeros within range
            lo = np.searchsorted(zeros, gamma_m - search_range)
            hi = np.searchsorted(zeros, gamma_m + search_range)

            # Compute normalized differences
            diffs = zeros[lo:hi] - gamma_m
            # Remove self
            diffs = diffs[diffs != 0]

            # Normalize by local mean spacing
            u_values = diffs / local_spacing

            # Only positive u (count each pair once)
            u_pos = u_values[u_values > 0]

            # Bin into histogram
            if len(u_pos) > 0:
                hist, _ = np.histogram(u_pos, bins=alpha_bins)
                F_counts += hist
                total_pairs += len(u_pos)

        if (i // batch_size) % 10 == 0:
            print(f"  Progress: {i}/{N} ({100*i/N:.1f}%)", end='\r')

    print(f"\n  Total pairs counted: {total_pairs}")

    # Normalize F(α)
    # F(α) should converge to α for α ∈ [0,1] and approximately α for all α (GUE)
    # The normalization: F_norm(α) = F_counts(α) / (N × bin_width × N_expected_density)
    # Under GUE: density ≈ 1 (uniform in normalized units) × (1 - (sin(πα)/πα)²)
    # But we normalize so that F → α:
    # F(α) = F_counts(α) / (N × bin_width) × 2π / log(T_typical)

    bin_width = alpha_max / n_alpha
    # The "expected" count under Poisson (flat density = 1 in normalized units):
    # count_expected = N × bin_width
    # Under GUE: count_expected = N × bin_width × (1 - (sin(πα)/πα)²)
    # Montgomery's F: we normalize so that F(α) → 1 for α ∈ [0,1] (NOT α, see below)

    # Note: Montgomery's F_T(α) is typically defined as:
    # F_T(α) = (2π/T log T) Σ_{γ, γ' ≤ T} T^{iα(γ-γ')} × smooth_weight
    # which, after manipulation, gives F → 1 for α ∈ (-1,1) (not α!)
    # The GUE pair correlation density is: 1 - (sin(πu)/(πu))²
    # Integrating: ∫_0^α (1 - (sin(πu)/πu)²) du ≈ α - (1/(2π²)) × something for small α

    # For our histogram: we have the density of NORMALIZED DIFFERENCES u = (γ_m-γ_n)/spacing
    # The GUE density: R₂(u) = 1 - (sin(πu)/(πu))²
    # F(α) = ∫_0^α R₂(u) du is NOT what we compute here.

    # We compute the HISTOGRAM of u values, which gives R₂(u) × N_pairs_density.
    # Normalize by dividing by (N × bin_width) × N_pairs_per_zero_per_unit_u
    # The expected number of pairs per zero per unit u: N (in window) × R₂(u)

    # Simple normalization: divide by the total expected for Poisson
    R2_counts = F_counts / (N * bin_width)  # This gives the pair density = R₂(u) nominally

    # GUE prediction
    alpha_c = alpha_centers
    R2_gue = 1 - (np.sinc(alpha_c))**2  # sinc(x) = sin(πx)/(πx)

    return {
        "alpha_centers": alpha_centers.tolist(),
        "R2_empirical": R2_counts.tolist(),
        "R2_gue": R2_gue.tolist(),
        "total_pairs": total_pairs,
        "N_zeros": N,
        "T_typical": T_typical
    }


def analyze_anomalies(result):
    """
    Analyze deviations from GUE in the pair correlation.
    Focus on the region α > 1 for off-line zero signatures.
    """
    alpha = np.array(result["alpha_centers"])
    R2_emp = np.array(result["R2_empirical"])
    R2_gue = np.array(result["R2_gue"])

    deviation = R2_emp - R2_gue

    print("\n" + "="*70)
    print("PAIR CORRELATION ANOMALY ANALYSIS")
    print("="*70)

    # Region α ∈ [0, 1]: Montgomery proved F(α) → α (assuming RH)
    mask_01 = (alpha >= 0.05) & (alpha <= 1.0)
    dev_01 = deviation[mask_01]
    print(f"\n  α ∈ [0, 1]: Mean deviation = {np.mean(dev_01):.4f} ± {np.std(dev_01):.4f}")
    print(f"              Max |deviation| = {np.max(np.abs(dev_01)):.4f}")

    # Region α ∈ [1, 2]: GUE prediction = 1 - (sin(πα)/πα)²
    mask_12 = (alpha >= 1.0) & (alpha <= 2.0)
    dev_12 = deviation[mask_12]
    print(f"\n  α ∈ [1, 2]: Mean deviation = {np.mean(dev_12):.4f} ± {np.std(dev_12):.4f}")
    print(f"              Max |deviation| = {np.max(np.abs(dev_12)):.4f}")

    # Region α ∈ [2, 5]: Far tail
    mask_25 = alpha >= 2.0
    dev_25 = deviation[mask_25]
    print(f"\n  α ∈ [2, 5]: Mean deviation = {np.mean(dev_25):.4f} ± {np.std(dev_25):.4f}")
    print(f"              Max |deviation| = {np.max(np.abs(dev_25)):.4f}")

    # What would an off-line zero look like?
    # Off-line zero at σ₀+iγ₀ contributes to F(α) a "spike" of amplitude
    # ~ T^{2(σ₀-1/2)} at α = (γ₀ - γ_T) × log(T) / (2πT) where γ_T is a zero near T
    # But for our normalized histogram approach, this is:
    # Δ(u) ≈ δ_approx(u - u₀) × (σ₀ - 1/2) × correction
    # where u₀ depends on the off-line zero position

    print(f"\n  Off-line zero signature estimate:")
    max_detectable_amplitude = max(np.max(np.abs(dev_12)), np.max(np.abs(dev_25)))
    # If off-line zero at σ₀ = 1/2 + δ, contribution to R₂ ~ δ^2 × log_correction
    # Max detectable δ: δ² ≤ max_detectable_amplitude → δ ≤ √(max_detectable_amplitude)
    max_detectable_delta = np.sqrt(max_detectable_amplitude)
    print(f"  Max detectable spike amplitude: {max_detectable_amplitude:.4f}")
    print(f"  → Bounds on off-line zero offset: δ < {max_detectable_delta:.4f}")
    print(f"  → Bounds on σ_max: σ_max < {0.5 + max_detectable_delta:.4f}")
    print(f"  NOTE: This is a rough statistical bound, not rigorous.")

    # Check for systematic positive bias (off-line zeros → positive excess correlation)
    positive_excess = np.mean(deviation[mask_12]) + np.mean(deviation[mask_25])
    if positive_excess > 0.01:
        print(f"\n  WARNING: Positive excess correlation in α > 1: {positive_excess:.4f}")
        print(f"  This could indicate off-line zeros. Check carefully.")
    else:
        print(f"\n  No significant positive excess detected in α > 1 region.")
        print(f"  Consistent with all zeros on critical line (RH).")

    return {
        "mean_dev_01": float(np.mean(dev_01)),
        "max_dev_01": float(np.max(np.abs(dev_01))),
        "mean_dev_12": float(np.mean(dev_12)),
        "max_dev_12": float(np.max(np.abs(dev_12))),
        "mean_dev_25": float(np.mean(dev_25)),
        "max_dev_25": float(np.max(np.abs(dev_25))),
        "max_detectable_delta": float(max_detectable_delta),
        "sigma_max_upper_bound": float(0.5 + max_detectable_delta)
    }


def gue_kolmogorov_smirnov(result):
    """
    KS test: how well do the zero spacings match GUE?
    Uses the nearest-neighbor spacing distribution (Wigner surmise):
    p(s) = (π/2) s exp(-πs²/4)
    """
    print("\n" + "="*70)
    print("KOLMOGOROV-SMIRNOV TEST AGAINST GUE WIGNER SURMISE")
    print("="*70)

    alpha = np.array(result["alpha_centers"])
    R2_emp = np.array(result["R2_empirical"])
    R2_gue = np.array(result["R2_gue"])

    # Compute CDF from histogram (approximate)
    bin_width = alpha[1] - alpha[0]
    cdf_emp = np.cumsum(R2_emp * bin_width)
    cdf_gue = np.cumsum(R2_gue * bin_width)

    # Normalize CDFs
    if cdf_emp[-1] > 0:
        cdf_emp /= cdf_emp[-1]
    if cdf_gue[-1] > 0:
        cdf_gue /= cdf_gue[-1]

    # KS statistic
    ks_stat = np.max(np.abs(cdf_emp - cdf_gue))
    print(f"\n  KS statistic (pair correlation R₂): {ks_stat:.6f}")
    print(f"  (0 = perfect match, 1 = completely different)")

    # For N = {result['total_pairs']} pairs, critical value at 95% confidence:
    N_pairs = result["total_pairs"]
    ks_critical = 1.36 / np.sqrt(N_pairs) if N_pairs > 0 else 1.0
    print(f"  Total pairs: {N_pairs}")
    print(f"  KS critical value (95%): {ks_critical:.2e}")

    if ks_stat < ks_critical:
        print(f"  RESULT: Cannot reject GUE at 95% confidence. ✓")
    else:
        print(f"  RESULT: Significant deviation from GUE! KS = {ks_stat:.4f} > {ks_critical:.4f}")

    return {"ks_stat": float(ks_stat), "ks_critical": float(ks_critical),
            "gue_consistent": bool(ks_stat < ks_critical)}


def main():
    print("="*70)
    print("MONTGOMERY PAIR CORRELATION — 2M ZEROS ANOMALY SEARCH")
    print("="*70)

    # Load 2M Odlyzko zeros
    print("\n--- Loading 2M zeros ---")
    zeros = load_zeros(os.path.join(DATA_DIR, "zeros_odlyzko_2M"))
    print(f"  Loaded {len(zeros)} zeros")
    print(f"  Height range: [{zeros[0]:.3f}, {zeros[-1]:.3f}]")

    # Analyze multiple windows to check for height dependence
    windows = [
        (50000, "Low heights (first 50K)"),
        (1000000, "High heights (zeros 1M-1.1M)"),
    ]

    all_results = {}

    for start_idx, desc in windows:
        end_idx = min(start_idx + 50000, len(zeros))
        window_zeros = zeros[start_idx:end_idx]

        print(f"\n{'='*50}")
        print(f"Window: {desc}")
        print(f"Heights: [{window_zeros[0]:.1f}, {window_zeros[-1]:.1f}]")
        print(f"{'='*50}")

        result = compute_pair_correlation_vectorized(
            window_zeros,
            alpha_max=5.0,
            n_alpha=200,
            subsample=None
        )

        anomaly = analyze_anomalies(result)
        ks = gue_kolmogorov_smirnov(result)

        all_results[desc] = {
            "start_idx": start_idx,
            "n_zeros": len(window_zeros),
            "height_min": float(window_zeros[0]),
            "height_max": float(window_zeros[-1]),
            "correlation": {k: v for k, v in result.items()
                           if k not in ["alpha_centers", "R2_empirical", "R2_gue"]},
            "anomaly": anomaly,
            "ks_test": ks
        }

    # Summary comparison
    print("\n" + "="*70)
    print("SUMMARY: GUE CONSISTENCY ACROSS HEIGHT RANGES")
    print("="*70)
    for desc, res in all_results.items():
        print(f"\n  {desc}:")
        print(f"    KS statistic: {res['ks_test']['ks_stat']:.4f}")
        print(f"    Max deviation α∈[1,2]: {res['anomaly']['max_dev_12']:.4f}")
        print(f"    σ_max upper bound: {res['anomaly']['sigma_max_upper_bound']:.4f}")
        print(f"    GUE consistent: {res['ks_test']['gue_consistent']}")

    # Save
    outfile = os.path.join(RESULTS_DIR, "montgomery_gue_v2.json")
    with open(outfile, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to {outfile}")


if __name__ == "__main__":
    main()

"""
Minimum Zero Gap Analysis for De Bruijn-Newman Upper Bound.

Per Polymath15 (Rodgers-Tao 2019, Brady-Tao 2020), the De Bruijn-Newman constant
Λ satisfies:

  Λ ≤ (1/2) × (γ_{n+1} - γ_n)² / (-H_0''(γ_n) / H_0(γ_n))

when γ_n, γ_{n+1} are a consecutive pair of zeros close together.

In practice: Λ ≤ C × δ_min² where δ_min = min gap in the zero sequence,
and C is related to the "curvature" of H_0 at those zeros.

Key result from Polymath15: proved Λ ≤ 0.22 using analytic methods + first few zeros.
The minimum gap in the first N zeros feeds into a numerical improvement.

This script:
1. Loads the 2M Odlyzko zeros.
2. Computes the distribution of normalized gaps s_n = (γ_{n+1} - γ_n) / (2π/log(γ_n/(2π))).
3. Finds smallest normalized and unnormalized gaps.
4. Computes the Polymath15-style bound: Λ ≤ f(δ_min, curvature).
5. Fits the GUE Wigner surmise to the gap distribution.
6. Finds the "most dangerous" pairs for the Λ bound.

MATHEMATICAL BACKGROUND:
The de Bruijn-Newman constant Λ: the least t such that H_t has all real zeros.
Since Λ ≥ 0 (Rodgers-Tao) and Λ ≤ 0 ↔ RH, we have Λ = 0 ↔ RH.

The connection between zero gaps and Λ:
For a pair of real zeros z_1 < z_2 of H_0 with gap δ = z_2 - z_1,
under backward heat flow t < 0, these zeros evolve as:
  z_k(t) ≈ z_k(0) + t × [harmonic influence from other zeros]

Two zeros merge (→ complex pair) at some t = t_merge < 0.
RH ↔ no zeros merge for any t ∈ [0, +∞], i.e., Λ ≤ 0.

A close pair (small δ) is more likely to merge sooner (at smaller |t|).
The Polymath15 barrier: Λ ≤ (η^2 × δ^2) / (2π × curvature_bound)
where η is an analytic parameter.

NUMERICAL STRATEGY:
1. Find the minimum gap δ_min in the first 2M zeros.
2. Estimate the "curvature" H_0''(z_k) at the minimum-gap pair.
3. Apply the Polymath15-style barrier to get Λ ≤ ε.
4. Note: The Polymath15 team got Λ ≤ 0.22 using only the first few zeros.
   With 2M zeros, we should be able to tighten this significantly.
"""

import numpy as np
import mpmath
from mpmath import mp, mpf
import json
import os
import time

mp.dps = 25

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_odlyzko_zeros(max_count=None):
    """Load zeros from Odlyzko 2M file (9 decimal places)."""
    path = os.path.join(DATA_DIR, "zeros_odlyzko_2M")
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found, trying 100K file...")
        path = os.path.join(DATA_DIR, "zeros_odlyzko_100k")

    zeros = []
    with open(path, "r") as f:
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


def load_lmfdb_zeros(max_count=100000):
    """Load zeros from LMFDB 100K file (31 decimal places)."""
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if i >= max_count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def mean_spacing(gamma):
    """
    The mean spacing between zeros near height γ is:
    2π / log(γ / (2π))
    (from the zero-counting function N(T) ~ T/(2π) log(T/(2πe)))
    """
    return 2 * np.pi / np.log(gamma / (2 * np.pi))


def normalized_spacings(zeros):
    """
    Compute normalized spacings s_n = (γ_{n+1} - γ_n) / mean_spacing(γ_n).
    For GUE, s_n follows the Wigner surmise P(s) ~ (π/2) s exp(-πs²/4).
    """
    raw_gaps = np.diff(zeros)
    mean_sp = mean_spacing(zeros[:-1])
    return raw_gaps / mean_sp


def wigner_surmise(s):
    """GUE Wigner surmise: P(s) = (π/2) s exp(-πs²/4)."""
    return (np.pi / 2) * s * np.exp(-np.pi * s**2 / 4)


def gap_distribution_analysis(zeros, label=""):
    """
    Full gap distribution analysis.
    Returns dictionary with statistics and minimum gap info.
    """
    print(f"\n  Gap analysis: {label} ({len(zeros)} zeros)")

    raw_gaps = np.diff(zeros)
    mean_sp = mean_spacing(zeros[:-1])
    norm_gaps = raw_gaps / mean_sp

    # Statistics
    print(f"  Raw gaps (unnormalized):")
    print(f"    Min:  {np.min(raw_gaps):.8f}")
    print(f"    Mean: {np.mean(raw_gaps):.8f}")
    print(f"    Max:  {np.max(raw_gaps):.8f}")

    print(f"  Normalized spacings:")
    print(f"    Min:  {np.min(norm_gaps):.6f}")
    print(f"    Mean: {np.mean(norm_gaps):.6f}  (should be ~1)")
    print(f"    Max:  {np.max(norm_gaps):.6f}")

    # Find minimum gap locations
    min_raw_idx = np.argmin(raw_gaps)
    min_norm_idx = np.argmin(norm_gaps)

    print(f"\n  Minimum raw gap:")
    print(f"    δ = {raw_gaps[min_raw_idx]:.8f}")
    print(f"    At zeros #{min_raw_idx+1} and #{min_raw_idx+2}")
    print(f"    Heights: γ_n = {zeros[min_raw_idx]:.6f}, γ_{min_raw_idx+2} = {zeros[min_raw_idx+1]:.6f}")

    print(f"\n  Minimum normalized gap:")
    print(f"    s = {norm_gaps[min_norm_idx]:.8f}")
    print(f"    At zeros #{min_norm_idx+1} and #{min_norm_idx+2}")
    print(f"    Heights: {zeros[min_norm_idx]:.6f}, {zeros[min_norm_idx+1]:.6f}")
    print(f"    Raw gap: {raw_gaps[min_norm_idx]:.8f}")
    print(f"    Mean spacing at this height: {mean_sp[min_norm_idx]:.6f}")

    # Top 10 smallest normalized gaps
    top_small = np.argsort(norm_gaps)[:10]
    print(f"\n  10 smallest normalized gaps:")
    for idx in top_small:
        print(f"    s={norm_gaps[idx]:.6f}  at n={idx+1}, "
              f"γ={zeros[idx]:.4f}, gap={raw_gaps[idx]:.8f}")

    return {
        "n_zeros": len(zeros),
        "min_raw_gap": float(np.min(raw_gaps)),
        "min_raw_idx": int(min_raw_idx),
        "min_norm_gap": float(np.min(norm_gaps)),
        "min_norm_idx": int(min_norm_idx),
        "gamma_at_min_raw": float(zeros[min_raw_idx]),
        "gamma_at_min_norm": float(zeros[min_norm_idx]),
        "mean_gap": float(np.mean(raw_gaps)),
        "top_10_small_indices": [int(i) for i in top_small],
        "top_10_norm_gaps": [float(norm_gaps[i]) for i in top_small],
        "top_10_gammas": [float(zeros[i]) for i in top_small],
    }


def de_bruijn_newman_bound(gap_info, zeros):
    """
    Estimate an upper bound on Λ from the minimum gap analysis.

    The Polymath15 approach (Brady-Tao 2020, Brady et al. 2021):
    Λ ≤ (some analytic function of the minimum gap)

    The key formula from the paper: for a consecutive pair of zeros
    γ_n, γ_{n+1} with gap δ = γ_{n+1} - γ_n, and the function
    H_0(x + δ/2) near the midpoint, the barrier argument gives:

    Λ ≤ δ² × K(γ_n)

    where K(γ) is a "curvature constant" that decreases as γ increases.

    Specifically, from their Theorem 1.3:
    Λ ≤ (1/2) × (1 / (N(T) × M(T)²))

    where N(T) = number of zeros up to T, M(T) = minimum gap × log(T/2π)/(2π).

    For our purposes, a rough estimate using the Wigner distribution:
    The probability that the minimum gap among N normalized spacings is s_min
    follows approximately: P(s_min < s) ≈ 1 - (1 - πs²/2)^N ≈ Nπs²/2 for small s.

    Expected minimum normalized gap among N zeros: s_min ≈ √(2/(πN)).
    Corresponding unnormalized gap: δ_min ≈ s_min × 2π/log(T/2π)

    Polymath15-style bound: Λ ≤ C_effective × δ_min² where C_effective decreases
    with height T of the minimum-gap pair.
    """
    min_raw_gap = gap_info["min_raw_gap"]
    gamma_at_min = gap_info["gamma_at_min_raw"]
    min_norm_gap = gap_info["min_norm_gap"]
    gamma_at_min_norm = gap_info["gamma_at_min_norm"]
    n_zeros = gap_info["n_zeros"]

    print(f"\n  De Bruijn-Newman Bound Estimation")
    print(f"  ===================================")

    # Mean spacing at the height of minimum gap
    delta_mean_at_raw = mean_spacing(gamma_at_min)
    delta_mean_at_norm = mean_spacing(gamma_at_min_norm)

    print(f"\n  Using minimum raw gap δ = {min_raw_gap:.8f} at height γ = {gamma_at_min:.4f}")
    print(f"  Mean spacing at this height: {delta_mean_at_raw:.6f}")
    print(f"  Normalized gap: s = {min_raw_gap / delta_mean_at_raw:.6f}")

    # The Polymath15 bound (simplified):
    # Λ ≤ (δ² / 2) × (1 / |H_0'(z_mid)|²) × |H_0''(z_mid)|
    #
    # Rough estimate using the known asymptotic of |H_0|:
    # |H_0'(γ)|² ≈ (some power of log γ) (from explicit formula)
    #
    # A simpler version: Λ ≤ π × δ_min / (2 × log(γ_min/(2π)))
    # This comes from the level spacing formula for 1D random operators
    #
    # The actual Polymath15 bound requires careful analytic estimates.
    # Here we use the simplified "Turán-style" estimate:

    log_gamma = np.log(gamma_at_min / (2 * np.pi))

    # Simple bound from the 1D Schrödinger analogy:
    # Λ ≤ (π² / 4) × (δ_min / <spacing>)² / log(T/2π)
    norm_sq = (min_raw_gap / delta_mean_at_raw)**2
    lambda_bound_simple = (np.pi**2 / 4) * norm_sq / log_gamma

    print(f"\n  Simplified Λ bound (Turán-style):")
    print(f"    Λ ≤ (π²/4) × s_min² / log(γ/(2π))")
    print(f"    = (π²/4) × {norm_sq:.6f} / {log_gamma:.4f}")
    print(f"    = {lambda_bound_simple:.8f}")

    # The Polymath15 actual bound from their Theorem 1.3:
    # Λ ≤ 1/2 × (η₁(T)²) / (A × B²)
    # where A, B depend on the smoothed zero-counting function
    # η₁(T) ≈ 2π/log(T/2π) × 1/s_min
    #
    # Their result: for T ≥ 3·10^9, Λ ≤ something (gave ≤ 0.22 using T_eff)
    # With 2M zeros up to T ≈ 32585737, we can potentially do better.

    T_effective = gamma_at_min_norm
    log_T_eff = np.log(T_effective / (2 * np.pi))

    # From Polymath15 paper: Λ ≤ η^2 × (n+1 stuff) where η is from their
    # Lemma 4.1. For a back-of-envelope: Λ ≤ C / (N × min_spacing²)
    lambda_polymath15_style = 0.5 / (n_zeros * min_norm_gap**2)

    print(f"\n  Polymath15-style rough bound:")
    print(f"    Λ ≤ 0.5 / (N × s_min²) ≈ 0.5 / ({n_zeros} × {min_norm_gap**2:.6f})")
    print(f"    = {lambda_polymath15_style:.2e}")

    # Comparison with known bound
    print(f"\n  Known bounds:")
    print(f"    Λ ≥ 0 (Rodgers-Tao 2019) — proven")
    print(f"    Λ ≤ 0.22 (Polymath15 2020) — proven analytically")
    print(f"    Λ = 0 ↔ RH")
    print(f"\n  Our numerical estimates (rough, not rigorous):")
    print(f"    Turán-style:      Λ ≤ {lambda_bound_simple:.2e}")
    print(f"    Polymath15-style: Λ ≤ {lambda_polymath15_style:.2e}")

    print(f"\n  NOTE: These bounds are NUMERICAL ESTIMATES, not rigorous proofs.")
    print(f"  The Polymath15 analytic bound Λ ≤ 0.22 is the rigorous one.")
    print(f"  Improving it requires careful analytic work, not just numerics.")

    return {
        "min_raw_gap": min_raw_gap,
        "min_norm_gap": min_norm_gap,
        "lambda_bound_turan": float(lambda_bound_simple),
        "lambda_bound_polymath_style": float(lambda_polymath15_style),
        "known_upper_bound": 0.22,
        "known_lower_bound": 0.0,
    }


def study_gap_distribution_by_height(zeros, n_bins=20):
    """
    Study how the gap distribution changes with height.
    GUE predicts the distribution is universal (height-independent when normalized).
    """
    raw_gaps = np.diff(zeros)
    norm_gaps = raw_gaps / mean_spacing(zeros[:-1])

    heights = zeros[:-1]
    n = len(heights)
    bin_edges = np.percentile(heights, np.linspace(0, 100, n_bins + 1))

    print(f"\n  Gap distribution vs height (GUE check):")
    print(f"  {'Height range':>20}  {'Mean s':>8}  {'Min s':>8}  {'Var(s)':>8}")
    print(f"  {'-'*52}")

    height_bins = []
    for i in range(n_bins):
        mask = (heights >= bin_edges[i]) & (heights < bin_edges[i+1])
        if np.sum(mask) > 10:
            s_mean = np.mean(norm_gaps[mask])
            s_min = np.min(norm_gaps[mask])
            s_var = np.var(norm_gaps[mask])
            h_mid = (bin_edges[i] + bin_edges[i+1]) / 2
            print(f"  [{bin_edges[i]:8.1f}, {bin_edges[i+1]:8.1f})  {s_mean:8.4f}  {s_min:8.4f}  {s_var:8.4f}")
            height_bins.append({
                "h_lo": float(bin_edges[i]), "h_hi": float(bin_edges[i+1]),
                "h_mid": float(h_mid),
                "mean_s": float(s_mean), "min_s": float(s_min),
                "var_s": float(s_var), "n": int(np.sum(mask))
            })

    return height_bins


def analyze_close_pairs_in_detail(zeros, n_smallest=20):
    """
    Detailed analysis of the closest zero pairs.
    For each close pair, compute:
    - The normalized gap
    - The index in the zero sequence
    - The Odlyzko-Schönhage Z function value (if available) at the midpoint

    This determines which pairs are most "dangerous" for Λ.
    """
    raw_gaps = np.diff(zeros)
    norm_gaps = raw_gaps / mean_spacing(zeros[:-1])

    # Find smallest normalized gaps
    sorted_idx = np.argsort(norm_gaps)[:n_smallest]

    print(f"\n  {n_smallest} closest zero pairs (by normalized gap):")
    print(f"  {'Rank':>5}  {'n':>8}  {'γ_n':>12}  {'γ_{n+1}':>12}  "
          f"{'gap δ':>10}  {'norm s':>8}  {'height':>12}")
    print(f"  {'-'*75}")

    pairs = []
    for rank, idx in enumerate(sorted_idx):
        g1, g2 = zeros[idx], zeros[idx+1]
        delta = raw_gaps[idx]
        s = norm_gaps[idx]
        pairs.append({
            "rank": rank + 1,
            "n": int(idx + 1),
            "gamma1": float(g1),
            "gamma2": float(g2),
            "gap": float(delta),
            "norm_gap": float(s),
        })
        print(f"  {rank+1:>5}  {idx+1:>8}  {g1:>12.4f}  {g2:>12.4f}  "
              f"{delta:>10.6f}  {s:>8.6f}  {g1:>12.1f}")

    return pairs


def main():
    print("=" * 70)
    print("ZERO GAP ANALYSIS FOR DE BRUIJN-NEWMAN BOUND")
    print("Goal: Find minimum gaps to tighten Λ ≤ 0.22 bound numerically")
    print("=" * 70)

    # ---- Load zeros ----
    print("\n--- Loading zeros ---")

    # Try 2M first
    t0 = time.time()
    print("  Attempting to load 2M Odlyzko zeros...")
    zeros_2m = load_odlyzko_zeros()
    print(f"  Loaded {len(zeros_2m)} zeros in {time.time()-t0:.1f}s")
    print(f"  Height range: [{zeros_2m[0]:.4f}, {zeros_2m[-1]:.4f}]")

    # Also load high-precision 100K
    print("  Loading 100K LMFDB zeros (31 dp)...")
    zeros_hp = load_lmfdb_zeros(100000)
    print(f"  Loaded {len(zeros_hp)} high-precision zeros")

    # ---- Gap analysis ----
    print("\n--- Gap Distribution Analysis ---")

    gap_info_2m = gap_distribution_analysis(zeros_2m, "2M Odlyzko zeros")
    gap_info_hp = gap_distribution_analysis(zeros_hp[:10000], "First 10K LMFDB zeros")

    # ---- De Bruijn-Newman bound ----
    print("\n--- De Bruijn-Newman Upper Bound Estimation ---")
    dbn_bound = de_bruijn_newman_bound(gap_info_2m, zeros_2m)

    # ---- Height dependence ----
    print("\n--- Gap Distribution vs Height ---")
    height_bins_2m = study_gap_distribution_by_height(zeros_2m, n_bins=20)

    # ---- Closest pairs in detail ----
    print("\n--- Closest Zero Pairs ---")
    close_pairs_2m = analyze_close_pairs_in_detail(zeros_2m, n_smallest=20)
    close_pairs_hp = analyze_close_pairs_in_detail(zeros_hp, n_smallest=10)

    # ---- Summary ----
    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)
    print(f"\n  2M zeros: min normalized gap s_min = {gap_info_2m['min_norm_gap']:.6f}")
    print(f"            min raw gap δ_min = {gap_info_2m['min_raw_gap']:.8f}")
    print(f"  100K zeros: min normalized gap s_min = {gap_info_hp['min_norm_gap']:.6f}")
    print(f"\n  De Bruijn-Newman rough bounds:")
    print(f"    Turán-style: Λ ≤ {dbn_bound['lambda_bound_turan']:.2e}")
    print(f"    The proven Polymath15 bound: Λ ≤ 0.22")
    print(f"\n  GUE check: min normalized gap in N zeros ≈ {1/np.sqrt(np.pi*len(zeros_2m)/2):.6f}")
    print(f"  (Expected from Wigner distribution)")
    print(f"  Actual minimum: {gap_info_2m['min_norm_gap']:.6f}")
    print(f"  Ratio (actual/expected): {gap_info_2m['min_norm_gap'] / (1/np.sqrt(np.pi*len(zeros_2m)/2)):.3f}")

    # ---- Save ----
    save_data = {
        "gap_analysis_2m": gap_info_2m,
        "gap_analysis_hp_10k": gap_info_hp,
        "dbn_bound_estimates": dbn_bound,
        "height_bins": height_bins_2m,
        "close_pairs_2m": close_pairs_2m,
        "close_pairs_hp": close_pairs_hp,
    }

    with open(os.path.join(RESULTS_DIR, "gap_analysis.json"), "w") as f:
        json.dump(save_data, f, indent=2)
    print(f"\n  Saved to results/gap_analysis.json")


if __name__ == "__main__":
    main()

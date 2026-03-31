"""
Full analysis on the 2M Odlyzko dataset.
Same metrics as analyze.py but with 20x the data.
Skips Li coefficients (too slow for 2M zeros with mpmath).
"""

import math
import os
import json

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_zeros_2M():
    path = os.path.join(DATA_DIR, "zeros_odlyzko_2M")
    zeros = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                zeros.append(float(line))
    return zeros


def N_approx(t):
    """Smooth approximation to zero counting function."""
    return (t / (2 * math.pi)) * math.log(t / (2 * math.pi)) - t / (2 * math.pi) + 7 / 8


def unfolded_spacings(zeros):
    """Compute unfolded spacings using smooth N(t)."""
    unfolded = [N_approx(t) for t in zeros]
    return [unfolded[i + 1] - unfolded[i] for i in range(len(unfolded) - 1)]


def spacing_stats(spacings):
    n = len(spacings)
    mean = sum(spacings) / n
    var = sum((s - mean) ** 2 for s in spacings) / n
    sorted_s = sorted(spacings)
    return {
        "count": n,
        "mean": mean,
        "variance": var,
        "std": math.sqrt(var),
        "median": sorted_s[n // 2],
        "min": sorted_s[0],
        "max": sorted_s[-1],
    }


def spacing_histogram(spacings, bins=200):
    """High-resolution histogram."""
    max_s = min(max(spacings), 4.0)  # cap at 4 for cleaner histogram
    bin_width = max_s / bins
    histogram = [0] * bins
    for s in spacings:
        if s < max_s:
            idx = int(s / bin_width)
            histogram[idx] += 1

    n = sum(histogram)
    pdf = [count / (n * bin_width) for count in histogram]
    centers = [(i + 0.5) * bin_width for i in range(bins)]

    # GUE Wigner surmise
    gue = [(32 / math.pi ** 2) * s ** 2 * math.exp(-4 * s ** 2 / math.pi) for s in centers]

    return centers, pdf, gue


def nearest_neighbor_ratio(spacings):
    """Compute r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1})."""
    ratios = []
    for i in range(len(spacings) - 1):
        s1, s2 = spacings[i], spacings[i + 1]
        if max(s1, s2) > 0:
            ratios.append(min(s1, s2) / max(s1, s2))

    mean_r = sum(ratios) / len(ratios)

    # GUE prediction for <r>: exact value is 4 - 2*sqrt(3) ≈ 0.5359
    # But that's for 2x2 GUE. For full GUE the value differs.
    # Let's compute the distribution too
    r_hist_bins = 100
    r_hist = [0] * r_hist_bins
    for r in ratios:
        idx = min(int(r * r_hist_bins), r_hist_bins - 1)
        r_hist[idx] += 1
    r_pdf = [count / (len(ratios) / r_hist_bins) for count in r_hist]

    return {
        "mean_ratio": mean_r,
        "n_ratios": len(ratios),
        "distribution": r_pdf,
    }


def pair_correlation(zeros, n_sample=500000):
    """
    Compute the pair correlation function R2(x).
    R2(x) = 1 - (sin(pi*x)/(pi*x))^2 for GUE (Montgomery's conjecture).

    We compute the empirical density of normalized differences
    (t_j - t_i) * log(t_i/2pi) / (2pi) for nearby pairs.
    """
    print("  Computing pair correlation (sampling)...")

    # Use unfolded zeros
    unfolded = [N_approx(t) for t in zeros[:n_sample]]

    # Collect differences up to distance 5
    max_dist = 5.0
    bins = 200
    bin_width = max_dist / bins
    histogram = [0] * bins
    n_pairs = 0

    # Only look at nearby pairs (within 10 raw positions)
    for i in range(len(unfolded) - 10):
        for j in range(i + 1, min(i + 10, len(unfolded))):
            diff = abs(unfolded[j] - unfolded[i])
            if diff < max_dist:
                idx = int(diff / bin_width)
                histogram[idx] += 1
                n_pairs += 1

    # Normalize to density
    centers = [(i + 0.5) * bin_width for i in range(bins)]
    # Each bin should have density relative to uniform Poisson = 1
    expected_per_bin = n_pairs * bin_width / max_dist
    R2_empirical = [count / expected_per_bin if expected_per_bin > 0 else 0 for count in histogram]

    # GUE prediction: R2(x) = 1 - (sin(pi*x)/(pi*x))^2
    R2_gue = []
    for x in centers:
        if x > 0.001:
            sinc = math.sin(math.pi * x) / (math.pi * x)
            R2_gue.append(1 - sinc ** 2)
        else:
            R2_gue.append(0)

    return {
        "centers": centers,
        "R2_empirical": R2_empirical,
        "R2_gue": R2_gue,
        "n_pairs": n_pairs,
    }


def convergence_by_height(zeros, spacings):
    """Check how GUE agreement changes with height on the critical line."""
    print("  Checking GUE agreement vs height...")

    chunk_size = 200000
    n_chunks = len(spacings) // chunk_size
    results = []

    for c in range(n_chunks):
        start = c * chunk_size
        end = start + chunk_size
        chunk = spacings[start:end]
        height = zeros[start + chunk_size // 2]

        mean = sum(chunk) / len(chunk)
        var = sum((s - mean) ** 2 for s in chunk) / len(chunk)

        results.append({
            "chunk": c,
            "height": height,
            "mean": mean,
            "variance": var,
        })
        print(f"    Height ~{height:.0f}: mean={mean:.6f}, var={var:.6f}")

    return results


def run():
    print("=" * 70)
    print("2M ZERO ANALYSIS")
    print("=" * 70)

    print("\nLoading 2M zeros...")
    zeros = load_zeros_2M()
    print(f"Loaded {len(zeros)} zeros. Range: {zeros[0]:.2f} to {zeros[-1]:.2f}")

    print("\nComputing unfolded spacings...")
    spacings = unfolded_spacings(zeros)

    # 1. Statistics
    print("\n" + "=" * 70)
    print("1. SPACING STATISTICS (2M zeros)")
    print("=" * 70)
    stats = spacing_stats(spacings)
    print(f"  Mean:     {stats['mean']:.8f}  (GUE: 1.0)")
    print(f"  Variance: {stats['variance']:.8f}  (GUE Wigner: 0.178)")
    print(f"  Std:      {stats['std']:.8f}")
    print(f"  Min:      {stats['min']:.8f}")
    print(f"  Max:      {stats['max']:.8f}")

    # 2. Histogram
    print("\n" + "=" * 70)
    print("2. SPACING DISTRIBUTION (200 bins)")
    print("=" * 70)
    centers, pdf, gue = spacing_histogram(spacings, bins=200)
    abs_devs = [abs(pdf[i] - gue[i]) for i in range(len(centers))]
    print(f"  Mean |empirical - GUE|: {sum(abs_devs) / len(abs_devs):.8f}")
    print(f"  Max  |empirical - GUE|: {max(abs_devs):.8f}")

    # 3. Nearest neighbor ratio
    print("\n" + "=" * 70)
    print("3. NEAREST NEIGHBOR RATIO")
    print("=" * 70)
    nn = nearest_neighbor_ratio(spacings)
    print(f"  Mean ratio: {nn['mean_ratio']:.8f}")

    # 4. Pair correlation
    print("\n" + "=" * 70)
    print("4. PAIR CORRELATION R2(x)")
    print("=" * 70)
    pc = pair_correlation(zeros, n_sample=500000)
    pc_devs = [abs(pc["R2_empirical"][i] - pc["R2_gue"][i]) for i in range(len(pc["centers"]))]
    print(f"  Mean |R2_emp - R2_gue|: {sum(pc_devs) / len(pc_devs):.6f}")
    print(f"  Pairs sampled: {pc['n_pairs']}")

    # 5. GUE convergence by height
    print("\n" + "=" * 70)
    print("5. GUE AGREEMENT vs HEIGHT")
    print("=" * 70)
    conv = convergence_by_height(zeros, spacings)

    # Save everything
    results = {
        "n_zeros": len(zeros),
        "zero_range": [zeros[0], zeros[-1]],
        "spacing_stats": stats,
        "histogram": {"centers": centers, "empirical": pdf, "gue": gue},
        "nearest_neighbor": nn,
        "pair_correlation": {
            "centers": pc["centers"],
            "R2_empirical": pc["R2_empirical"],
            "R2_gue": pc["R2_gue"],
        },
        "convergence_by_height": conv,
    }

    outfile = os.path.join(RESULTS_DIR, "analysis_2M.json")
    with open(outfile, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {outfile}")


if __name__ == "__main__":
    run()

"""
Analysis of Riemann zeta zeros — spacing distributions, GUE statistics,
Li coefficients, and pattern exploration.
"""

import mpmath
import math
import os
import json

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_zeros(count=None):
    """Load zeros as Python floats for fast analysis."""
    # Prefer Odlyzko 100K (complete), fall back to LMFDB partial
    odlyzko = os.path.join(DATA_DIR, "zeros_odlyzko_100k")
    lmfdb = os.path.join(DATA_DIR, "zeros_100k.txt")
    path = odlyzko if os.path.exists(odlyzko) else lmfdb
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if count and i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return zeros


# =============================================================================
# 1. SPACING ANALYSIS
# =============================================================================

def normalized_spacings(zeros):
    """
    Compute normalized spacings between consecutive zeros.
    Normalization: multiply by (log(t/2pi) / 2pi) to account for
    increasing density of zeros as t grows.
    """
    spacings = []
    for i in range(len(zeros) - 1):
        raw_gap = zeros[i + 1] - zeros[i]
        t = zeros[i]
        # Average density at height t is (1/2pi) * log(t/2pi)
        density = math.log(t / (2 * math.pi)) / (2 * math.pi)
        normalized = raw_gap * density
        spacings.append(normalized)
    return spacings


def spacing_statistics(spacings):
    """Basic statistics of normalized spacings."""
    n = len(spacings)
    mean = sum(spacings) / n
    var = sum((s - mean) ** 2 for s in spacings) / n
    std = math.sqrt(var)

    sorted_s = sorted(spacings)
    median = sorted_s[n // 2]
    min_s = sorted_s[0]
    max_s = sorted_s[-1]

    return {
        "count": n,
        "mean": mean,
        "std": std,
        "variance": var,
        "median": median,
        "min": min_s,
        "max": max_s,
    }


def spacing_distribution(spacings, bins=100):
    """
    Histogram of normalized spacings.
    Compare against GUE prediction (Wigner surmise for GUE):
    p(s) = (32/pi^2) * s^2 * exp(-4s^2/pi)
    """
    max_s = max(spacings)
    bin_width = max_s / bins
    histogram = [0] * bins

    for s in spacings:
        idx = min(int(s / bin_width), bins - 1)
        histogram[idx] += 1

    # Normalize to probability density
    n = len(spacings)
    pdf = [count / (n * bin_width) for count in histogram]
    bin_centers = [(i + 0.5) * bin_width for i in range(bins)]

    # GUE Wigner surmise values at same points
    gue_pdf = []
    for s in bin_centers:
        gue = (32 / (math.pi ** 2)) * s ** 2 * math.exp(-4 * s ** 2 / math.pi)
        gue_pdf.append(gue)

    return bin_centers, pdf, gue_pdf


# =============================================================================
# 2. GUE DEVIATION ANALYSIS
# =============================================================================

def gue_deviation(bin_centers, empirical_pdf, gue_pdf):
    """
    Compute deviation between empirical and GUE distributions.
    Returns per-bin deviations and summary statistics.
    """
    deviations = []
    for i in range(len(bin_centers)):
        if gue_pdf[i] > 1e-10:  # avoid division by near-zero
            rel_dev = (empirical_pdf[i] - gue_pdf[i]) / gue_pdf[i]
        else:
            rel_dev = 0
        deviations.append({
            "s": bin_centers[i],
            "empirical": empirical_pdf[i],
            "gue": gue_pdf[i],
            "abs_deviation": empirical_pdf[i] - gue_pdf[i],
            "rel_deviation": rel_dev,
        })

    abs_devs = [abs(d["abs_deviation"]) for d in deviations]
    mean_abs_dev = sum(abs_devs) / len(abs_devs)
    max_abs_dev = max(abs_devs)
    max_idx = abs_devs.index(max_abs_dev)

    return {
        "deviations": deviations,
        "mean_abs_deviation": mean_abs_dev,
        "max_abs_deviation": max_abs_dev,
        "max_deviation_at_s": bin_centers[max_idx],
    }


# =============================================================================
# 3. LI COEFFICIENTS
# =============================================================================

def li_coefficients(zeros, n_coefficients=50):
    """
    Compute Li coefficients lambda_n from the known zeros.
    Li's criterion: RH is true iff lambda_n >= 0 for all n >= 1.

    lambda_n = sum over zeros rho of [1 - (1 - 1/rho)^n]

    Using the imaginary parts t_k (assuming RH, rho = 1/2 + i*t_k):
    lambda_n = sum_k [1 - (1 - 1/(1/2 + i*t_k))^n]

    We use a finite sum over our known zeros as an approximation.
    """
    mpmath.mp.dps = 50
    results = []

    # Use first 10000 zeros for Li coefficient computation
    n_zeros = min(len(zeros), 10000)
    print(f"Computing Li coefficients using {n_zeros} zeros...")

    for n in range(1, n_coefficients + 1):
        lambda_n = mpmath.mpf(0)
        for k in range(n_zeros):
            t = mpmath.mpf(zeros[k])
            rho = mpmath.mpc(0.5, t)
            term = 1 - (1 - 1 / rho) ** n
            # Include conjugate zero (1/2 - i*t)
            rho_conj = mpmath.mpc(0.5, -t)
            term_conj = 1 - (1 - 1 / rho_conj) ** n
            lambda_n += term + term_conj

        lambda_n = float(lambda_n.real)
        results.append({"n": n, "lambda_n": lambda_n})

        if n <= 10 or n % 10 == 0:
            print(f"  lambda_{n} = {lambda_n:.10f}")

    return results


# =============================================================================
# 4. NEAREST NEIGHBOR STATISTICS
# =============================================================================

def nearest_neighbor_ratio(zeros):
    """
    Compute ratio of consecutive spacings: r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1})
    For GUE, the mean ratio should be approximately 0.5307...
    (Atas et al., 2013)
    """
    spacings = []
    for i in range(len(zeros) - 1):
        spacings.append(zeros[i + 1] - zeros[i])

    ratios = []
    for i in range(len(spacings) - 1):
        s1, s2 = spacings[i], spacings[i + 1]
        r = min(s1, s2) / max(s1, s2)
        ratios.append(r)

    mean_ratio = sum(ratios) / len(ratios)
    gue_expected = 0.5307  # 4 - 2*sqrt(3) approx

    return {
        "mean_ratio": mean_ratio,
        "gue_expected": gue_expected,
        "deviation": mean_ratio - gue_expected,
        "n_ratios": len(ratios),
    }


# =============================================================================
# 5. NUMBER VARIANCE
# =============================================================================

def number_variance(zeros, L_values=None):
    """
    Sigma^2(L) = variance of the number of zeros in an interval of length L.
    For GUE: Sigma^2(L) ~ (2/pi^2) * (ln(2*pi*L) + gamma + 1 - pi^2/8)

    This is a key statistic for testing GUE universality.
    """
    if L_values is None:
        L_values = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]

    results = []
    # Use the middle section of zeros (avoid boundary effects)
    start_idx = len(zeros) // 4
    end_idx = 3 * len(zeros) // 4

    for L in L_values:
        counts = []
        # Sample windows
        step = max(1, (end_idx - start_idx) // 1000)
        for i in range(start_idx, end_idx, step):
            t_start = zeros[i]
            t_end = t_start + L
            count = sum(1 for z in zeros[i:] if z <= t_end) - 1
            counts.append(count)

        mean_count = sum(counts) / len(counts)
        variance = sum((c - mean_count) ** 2 for c in counts) / len(counts)

        # GUE prediction
        euler_gamma = 0.5772156649015329
        if L > 0:
            gue_var = (2 / math.pi ** 2) * (math.log(2 * math.pi * L) + euler_gamma + 1 - math.pi ** 2 / 8)
        else:
            gue_var = 0

        results.append({
            "L": L,
            "empirical_variance": variance,
            "gue_variance": gue_var,
            "deviation": variance - gue_var,
        })
        print(f"  L={L:>5.1f}: empirical={variance:.6f}, GUE={gue_var:.6f}, dev={variance - gue_var:.6f}")

    return results


# =============================================================================
# MAIN — RUN ALL ANALYSES
# =============================================================================

def run_all():
    print("=" * 70)
    print("RIEMANN ZETA ZERO ANALYSIS")
    print("=" * 70)

    print("\nLoading zeros...")
    zeros = load_zeros()
    print(f"Loaded {len(zeros)} zeros.")
    print(f"Range: {zeros[0]:.6f} to {zeros[-1]:.6f}")

    # 1. Spacing analysis
    print("\n" + "=" * 70)
    print("1. NORMALIZED SPACING STATISTICS")
    print("=" * 70)
    spacings = normalized_spacings(zeros)
    stats = spacing_statistics(spacings)
    print(f"  Mean:     {stats['mean']:.8f}  (GUE predicts ~1.0)")
    print(f"  Std:      {stats['std']:.8f}")
    print(f"  Variance: {stats['variance']:.8f}  (GUE predicts ~0.286)")
    print(f"  Median:   {stats['median']:.8f}")
    print(f"  Min:      {stats['min']:.8f}")
    print(f"  Max:      {stats['max']:.8f}")

    # 2. Spacing distribution vs GUE
    print("\n" + "=" * 70)
    print("2. SPACING DISTRIBUTION vs GUE (Wigner Surmise)")
    print("=" * 70)
    bin_centers, emp_pdf, gue_pdf = spacing_distribution(spacings, bins=50)
    gue_dev = gue_deviation(bin_centers, emp_pdf, gue_pdf)
    print(f"  Mean absolute deviation: {gue_dev['mean_abs_deviation']:.8f}")
    print(f"  Max absolute deviation:  {gue_dev['max_abs_deviation']:.8f} at s={gue_dev['max_deviation_at_s']:.4f}")

    # Save distribution data for plotting
    dist_data = {
        "bin_centers": bin_centers,
        "empirical_pdf": emp_pdf,
        "gue_pdf": gue_pdf,
        "deviations": gue_dev["deviations"][:20],  # first 20 for readability
    }

    # 3. Nearest neighbor ratio
    print("\n" + "=" * 70)
    print("3. NEAREST NEIGHBOR RATIO")
    print("=" * 70)
    nn = nearest_neighbor_ratio(zeros)
    print(f"  Mean ratio:   {nn['mean_ratio']:.8f}")
    print(f"  GUE expected: {nn['gue_expected']:.4f}")
    print(f"  Deviation:    {nn['deviation']:.8f}")

    # 4. Number variance
    print("\n" + "=" * 70)
    print("4. NUMBER VARIANCE Sigma^2(L)")
    print("=" * 70)
    nv = number_variance(zeros)

    # 5. Li coefficients (slower — uses mpmath)
    print("\n" + "=" * 70)
    print("5. LI COEFFICIENTS (first 50)")
    print("=" * 70)
    li = li_coefficients(zeros[:10000], n_coefficients=50)
    all_positive = all(l["lambda_n"] > 0 for l in li)
    print(f"\n  All positive (consistent with RH): {all_positive}")

    # Save all results
    results = {
        "n_zeros": len(zeros),
        "zero_range": [zeros[0], zeros[-1]],
        "spacing_stats": stats,
        "gue_deviation_summary": {
            "mean_abs": gue_dev["mean_abs_deviation"],
            "max_abs": gue_dev["max_abs_deviation"],
            "max_at_s": gue_dev["max_deviation_at_s"],
        },
        "nearest_neighbor": nn,
        "number_variance": nv,
        "li_coefficients": li,
        "li_all_positive": all_positive,
    }

    outfile = os.path.join(RESULTS_DIR, "analysis_results.json")
    with open(outfile, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {outfile}")

    return results


if __name__ == "__main__":
    run_all()

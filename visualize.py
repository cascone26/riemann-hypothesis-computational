"""
Visualizations of Riemann zeta zero analysis.
Generates PNG plots for: spacing distribution, Li coefficients,
pair correlation, and zero distribution on the critical line.
"""

import matplotlib
matplotlib.use('Agg')  # headless
import matplotlib.pyplot as plt
import math
import json
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
PLOTS_DIR = os.path.join(os.path.dirname(__file__), "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)


def load_zeros(source="odlyzko_100k"):
    if source == "odlyzko_100k":
        path = os.path.join(DATA_DIR, "zeros_odlyzko_100k")
    else:
        path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                zeros.append(float(line))
    return zeros


def N_approx(t):
    return (t / (2 * math.pi)) * math.log(t / (2 * math.pi)) - t / (2 * math.pi) + 7 / 8


def plot_spacing_distribution(zeros):
    """Plot empirical spacing distribution vs GUE Wigner surmise."""
    print("Plotting spacing distribution...")

    unfolded = [N_approx(t) for t in zeros]
    spacings = [unfolded[i + 1] - unfolded[i] for i in range(len(unfolded) - 1)]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Histogram vs GUE
    ax1.hist(spacings, bins=150, density=True, alpha=0.7, color='steelblue',
             label='Empirical (100K zeros)', range=(0, 3.5))

    s_vals = [i * 0.01 for i in range(1, 350)]
    gue = [(32 / math.pi ** 2) * s ** 2 * math.exp(-4 * s ** 2 / math.pi) for s in s_vals]
    ax1.plot(s_vals, gue, 'r-', linewidth=2, label='GUE Wigner Surmise')

    ax1.set_xlabel('Normalized Spacing s', fontsize=12)
    ax1.set_ylabel('Probability Density', fontsize=12)
    ax1.set_title('Riemann Zeta Zero Spacings vs GUE Prediction', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.set_xlim(0, 3.5)

    # Residuals
    bin_edges = [i * 0.025 for i in range(140)]
    bin_centers = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]
    hist_counts = [0] * (len(bin_edges) - 1)
    bin_width = bin_edges[1] - bin_edges[0]

    for s in spacings:
        idx = int(s / bin_width)
        if 0 <= idx < len(hist_counts):
            hist_counts[idx] += 1

    n = len(spacings)
    emp_pdf = [c / (n * bin_width) for c in hist_counts]
    gue_at_centers = [(32 / math.pi ** 2) * s ** 2 * math.exp(-4 * s ** 2 / math.pi) for s in bin_centers]
    residuals = [emp_pdf[i] - gue_at_centers[i] for i in range(len(bin_centers))]

    ax2.bar(bin_centers, residuals, width=bin_width * 0.9, color='steelblue', alpha=0.7)
    ax2.axhline(y=0, color='red', linewidth=1)
    ax2.set_xlabel('Normalized Spacing s', fontsize=12)
    ax2.set_ylabel('Empirical - GUE', fontsize=12)
    ax2.set_title('Residuals (Empirical minus GUE Wigner Surmise)', fontsize=14)
    ax2.set_xlim(0, 3.5)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'spacing_distribution.png'), dpi=150)
    plt.close()
    print("  Saved spacing_distribution.png")


def plot_li_coefficients():
    """Plot Li coefficients and compare to asymptotic formula."""
    print("Plotting Li coefficients...")

    # Load from results
    path = os.path.join(RESULTS_DIR, "analysis_results.json")
    with open(path) as f:
        data = json.load(f)

    li = data["li_coefficients"]
    ns = [l["n"] for l in li]
    lambdas = [l["lambda_n"] for l in li]

    gamma = 0.5772156649015329
    asymptotic = []
    for n in ns:
        if n >= 2:
            a = (n / 2) * (math.log(n) + gamma - 1 - math.log(2 * math.pi))
        else:
            a = 0.023
        asymptotic.append(a)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    ax1.plot(ns, lambdas, 'b.-', markersize=4, label='Computed lambda_n')
    ax1.plot(ns, asymptotic, 'r--', linewidth=1.5, label='Keiper-Li asymptotic')
    ax1.set_xlabel('n', fontsize=12)
    ax1.set_ylabel('lambda_n', fontsize=12)
    ax1.set_title('Li Coefficients (all positive iff RH is true)', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.axhline(y=0, color='gray', linestyle=':', alpha=0.5)

    # Ratio
    ratios = [lambdas[i] / asymptotic[i] if asymptotic[i] != 0 and asymptotic[i] > 0.01 else 1
              for i in range(len(ns))]
    ax2.plot(ns[5:], ratios[5:], 'g.-', markersize=4)
    ax2.axhline(y=1.0, color='red', linestyle='--', linewidth=1)
    ax2.set_xlabel('n', fontsize=12)
    ax2.set_ylabel('lambda_n / asymptotic', fontsize=12)
    ax2.set_title('Ratio to Keiper-Li Asymptotic (should approach 1)', fontsize=14)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'li_coefficients.png'), dpi=150)
    plt.close()
    print("  Saved li_coefficients.png")


def plot_pair_correlation(zeros):
    """Plot the pair correlation function R2(x) vs GUE prediction."""
    print("Plotting pair correlation...")

    n_sample = min(len(zeros), 100000)
    unfolded = [N_approx(t) for t in zeros[:n_sample]]

    max_dist = 5.0
    bins = 200
    bin_width = max_dist / bins
    histogram = [0] * bins
    n_pairs = 0

    for i in range(len(unfolded) - 10):
        for j in range(i + 1, min(i + 10, len(unfolded))):
            diff = abs(unfolded[j] - unfolded[i])
            if diff < max_dist:
                idx = int(diff / bin_width)
                if idx < bins:
                    histogram[idx] += 1
                    n_pairs += 1

    centers = [(i + 0.5) * bin_width for i in range(bins)]
    expected_per_bin = n_pairs * bin_width / max_dist
    R2_emp = [count / expected_per_bin if expected_per_bin > 0 else 0 for count in histogram]

    R2_gue = []
    for x in centers:
        if x > 0.001:
            sinc = math.sin(math.pi * x) / (math.pi * x)
            R2_gue.append(1 - sinc ** 2)
        else:
            R2_gue.append(0)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(centers, R2_emp, 'b-', alpha=0.7, linewidth=1, label='Empirical R2(x)')
    ax.plot(centers, R2_gue, 'r--', linewidth=2, label='GUE: 1 - (sin(pi*x)/(pi*x))^2')
    ax.set_xlabel('x (unfolded spacing)', fontsize=12)
    ax.set_ylabel('R2(x)', fontsize=12)
    ax.set_title('Pair Correlation Function — Montgomery-Odlyzko Law', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0, 5)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'pair_correlation.png'), dpi=150)
    plt.close()
    print("  Saved pair_correlation.png")


def plot_zeros_on_critical_line(zeros):
    """Visualize the first ~200 zeros on the critical line."""
    print("Plotting zeros on critical line...")

    n_show = 200
    subset = zeros[:n_show]

    fig, ax = plt.subplots(figsize=(4, 16))
    for t in subset:
        ax.plot(0.5, t, 'b.', markersize=3)

    ax.set_xlim(0, 1)
    ax.set_xlabel('Re(s)', fontsize=12)
    ax.set_ylabel('Im(s)', fontsize=12)
    ax.set_title(f'First {n_show} Non-Trivial Zeros\n(all on Re(s) = 1/2)', fontsize=13)
    ax.axvline(x=0.5, color='red', linestyle='--', alpha=0.3, label='Critical line')
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'critical_line.png'), dpi=150)
    plt.close()
    print("  Saved critical_line.png")


def plot_spacing_convergence(zeros):
    """Show how variance converges to GUE with more zeros."""
    print("Plotting convergence...")

    sizes = [100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000]
    variances = []

    for n in sizes:
        subset = zeros[:n]
        unfolded = [N_approx(t) for t in subset]
        spacings = [unfolded[i + 1] - unfolded[i] for i in range(len(unfolded) - 1)]
        mean = sum(spacings) / len(spacings)
        var = sum((s - mean) ** 2 for s in spacings) / len(spacings)
        variances.append(var)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogx(sizes, variances, 'bo-', markersize=6, label='Empirical variance')
    ax.axhline(y=0.178, color='red', linestyle='--', linewidth=1.5,
               label='GUE Wigner surmise (0.178)')
    ax.set_xlabel('Number of zeros', fontsize=12)
    ax.set_ylabel('Spacing variance', fontsize=12)
    ax.set_title('Convergence of Spacing Variance to GUE', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'convergence.png'), dpi=150)
    plt.close()
    print("  Saved convergence.png")


def main():
    print("=" * 50)
    print("GENERATING VISUALIZATIONS")
    print("=" * 50)

    zeros = load_zeros("odlyzko_100k")
    print(f"Loaded {len(zeros)} zeros\n")

    plot_spacing_distribution(zeros)
    plot_li_coefficients()
    plot_pair_correlation(zeros)
    plot_zeros_on_critical_line(zeros)
    plot_spacing_convergence(zeros)

    print(f"\nAll plots saved to {PLOTS_DIR}/")


if __name__ == "__main__":
    main()

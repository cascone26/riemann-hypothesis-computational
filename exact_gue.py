"""
Exact GUE nearest-neighbor spacing distribution via Fredholm determinant
of the sine kernel, compared against empirical zeta zero spacings and
the Wigner surmise approximation.

Method (Gaudin-Mehta):
  E(s) = det(I - K_s)  where K_s is the sine kernel on [0, s]
  p(s) = d^2/ds^2 E(s)

The sine kernel: K(x,y) = sin(pi(x-y)) / (pi(x-y))
Discretized via Gauss-Legendre quadrature on [0, s].
"""

import numpy as np
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE, "data")
PLOTS_DIR = os.path.join(BASE, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Exact GUE spacing via Fredholm determinant of sine kernel
# ---------------------------------------------------------------------------

def sine_kernel_matrix(s, n_quad=128):
    """
    Build the n_quad x n_quad discretized sine kernel matrix on [0, s]
    using Gauss-Legendre quadrature.

    Returns det(I - K) where K_{ij} = sqrt(w_i) * K(x_i, x_j) * sqrt(w_j)
    with K(x,y) = sin(pi(x-y)) / (pi(x-y)).
    """
    # Gauss-Legendre on [-1, 1], then shift to [0, s]
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    # Transform to [0, s]
    x = 0.5 * s * (nodes + 1.0)
    w = 0.5 * s * weights

    # Build kernel matrix
    # K(x_i, x_j) = sin(pi(x_i - x_j)) / (pi(x_i - x_j))
    diff = x[:, None] - x[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        arg = np.pi * diff
        K = np.where(np.abs(diff) < 1e-14, 1.0, np.sin(arg) / arg)

    # Incorporate quadrature weights: symmetric form
    sw = np.sqrt(w)
    K_weighted = sw[:, None] * K * sw[None, :]

    return K_weighted


def fredholm_det(s, n_quad=128):
    """Compute E(s) = det(I - K_s) via eigenvalues of discretized kernel."""
    if s < 1e-14:
        return 1.0
    K = sine_kernel_matrix(s, n_quad)
    eigenvalues = np.linalg.eigvalsh(K)
    # det(I - K) = prod(1 - lambda_i)
    log_det = np.sum(np.log(np.maximum(1.0 - eigenvalues, 1e-300)))
    return np.exp(log_det)


def compute_exact_gue(s_values, n_quad=128):
    """
    Compute E(s), then obtain p(s) = d^2E/ds^2 via finite differences.
    Returns (E_array, pdf_array).
    """
    E = np.array([fredholm_det(s, n_quad) for s in s_values])
    return E


def exact_gue_pdf_from_E(s_values, E_values):
    """
    p(s) = d^2 E(s) / ds^2 computed via central finite differences.
    """
    ds = s_values[1] - s_values[0]
    # Second derivative via central differences
    pdf = np.zeros_like(E_values)
    pdf[1:-1] = (E_values[2:] - 2.0 * E_values[1:-1] + E_values[:-2]) / (ds * ds)
    # Forward/backward at endpoints
    if len(s_values) >= 3:
        pdf[0] = (E_values[2] - 2.0 * E_values[1] + E_values[0]) / (ds * ds)
        pdf[-1] = (E_values[-1] - 2.0 * E_values[-2] + E_values[-3]) / (ds * ds)
    return pdf


# ---------------------------------------------------------------------------
# 2. Wigner surmise (GUE approximation)
# ---------------------------------------------------------------------------

def wigner_surmise(s):
    """p_W(s) = (32/pi^2) * s^2 * exp(-4s^2/pi)"""
    return (32.0 / (np.pi**2)) * s**2 * np.exp(-4.0 * s**2 / np.pi)


# ---------------------------------------------------------------------------
# 3. Load zeros and compute unfolded spacings
# ---------------------------------------------------------------------------

def load_zeros():
    """Load Odlyzko 100K zeros."""
    path = os.path.join(DATA_DIR, "zeros_odlyzko_100k")
    zeros = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def unfold_zeros(zeros):
    """
    Unfold using smooth counting function:
    N(t) = (t/(2*pi)) * log(t/(2*pi)) - t/(2*pi) + 7/8
    Returns unfolded spacings (should have mean ~1).
    """
    t = zeros
    twopi = 2.0 * np.pi
    N = (t / twopi) * np.log(t / twopi) - t / twopi + 7.0 / 8.0
    spacings = np.diff(N)
    return spacings


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("EXACT GUE SPACING DISTRIBUTION vs ZETA ZEROS")
    print("=" * 60)

    # --- Compute exact GUE ---
    n_points = 500
    s_max = 4.0
    # Use a finer grid for accurate second derivative
    n_fine = 2000
    s_fine = np.linspace(0, s_max, n_fine)

    print(f"\nComputing Fredholm determinant E(s) at {n_fine} points...")
    print("(Gauss-Legendre quadrature with 128 nodes)")
    E_fine = compute_exact_gue(s_fine, n_quad=128)

    print("Computing exact GUE PDF via d^2E/ds^2...")
    pdf_fine = exact_gue_pdf_from_E(s_fine, E_fine)

    # Interpolate to the 500-point grid
    s_grid = np.linspace(0, s_max, n_points)
    E_grid = np.interp(s_grid, s_fine, E_fine)
    pdf_exact = np.interp(s_grid, s_fine, pdf_fine)

    # Compute CDF from E(s): CDF(s) = 1 - E(s) + integral of p from 0 to s
    # Actually, the gap probability E(s) = P(no eigenvalue in [0,s])
    # The nearest-neighbor CDF is: F(s) = 1 - E(s) ... wait, not exactly.
    # For GUE: F(s) = 1 - E(s) is the probability that the nearest spacing <= s
    # when p(s) = d^2E/ds^2. Actually p(s) = d^2E/ds^2 and integral of p = 1.
    # The CDF is just the cumulative integral of pdf.
    ds = s_grid[1] - s_grid[0]
    cdf_exact = np.cumsum(pdf_exact) * ds

    print(f"  PDF integral check: {np.trapezoid(pdf_exact, s_grid):.6f} (should be ~1)")
    print(f"  PDF max: {np.max(pdf_exact):.6f} at s={s_grid[np.argmax(pdf_exact)]:.4f}")

    # --- Exact GUE variance ---
    mean_exact = np.trapezoid(s_grid * pdf_exact, s_grid)
    var_exact = np.trapezoid(s_grid**2 * pdf_exact, s_grid) - mean_exact**2
    print(f"  Exact GUE mean: {mean_exact:.6f}")
    print(f"  Exact GUE variance: {var_exact:.6f}")

    # --- Wigner surmise ---
    pdf_wigner = wigner_surmise(s_grid)

    # --- Load and process zeros ---
    print("\nLoading 100K Odlyzko zeros...")
    zeros = load_zeros()
    print(f"  Loaded {len(zeros)} zeros, range [{zeros[0]:.3f}, {zeros[-1]:.3f}]")

    spacings = unfold_zeros(zeros)
    print(f"  {len(spacings)} unfolded spacings")
    print(f"  Mean spacing: {np.mean(spacings):.6f} (should be ~1)")
    print(f"  Variance: {np.var(spacings):.6f}")

    # --- Empirical histogram (density) ---
    n_bins = 150
    hist_counts, bin_edges = np.histogram(spacings, bins=n_bins, range=(0, s_max), density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # --- Compute deviations ---
    # Interpolate exact and Wigner onto bin centers
    pdf_exact_at_bins = np.interp(bin_centers, s_grid, pdf_exact)
    pdf_wigner_at_bins = wigner_surmise(bin_centers)

    mad_wigner = np.mean(np.abs(hist_counts - pdf_wigner_at_bins))
    mad_exact = np.mean(np.abs(hist_counts - pdf_exact_at_bins))

    print(f"\n--- Deviations ---")
    print(f"  Mean absolute deviation from Wigner surmise: {mad_wigner:.6f}")
    print(f"  Mean absolute deviation from exact GUE:      {mad_exact:.6f}")
    print(f"  Improvement factor: {mad_wigner / mad_exact:.2f}x")

    # --- Plot 1: Main comparison ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1]})

    ax1.bar(bin_centers, hist_counts, width=bin_edges[1]-bin_edges[0],
            alpha=0.4, color='steelblue', label='Zeta zeros (100K Odlyzko)')
    ax1.plot(s_grid, pdf_wigner, 'r--', linewidth=2, label='Wigner surmise (approx)')
    ax1.plot(s_grid, pdf_exact, 'k-', linewidth=2, label='Exact GUE (Fredholm det)')
    ax1.set_xlabel('Normalized spacing s', fontsize=13)
    ax1.set_ylabel('p(s)', fontsize=13)
    ax1.set_title('Zeta Zero Spacings vs GUE: Wigner Surmise and Exact Distribution', fontsize=14)
    ax1.legend(fontsize=12)
    ax1.set_xlim(0, s_max)
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)

    # Inset: zoom on peak region
    ax_inset = ax1.inset_axes([0.55, 0.5, 0.4, 0.4])
    peak_mask = (s_grid > 0.6) & (s_grid < 1.4)
    bin_mask = (bin_centers > 0.6) & (bin_centers < 1.4)
    ax_inset.bar(bin_centers[bin_mask], hist_counts[bin_mask],
                 width=bin_edges[1]-bin_edges[0], alpha=0.4, color='steelblue')
    ax_inset.plot(s_grid[peak_mask], pdf_wigner[peak_mask], 'r--', linewidth=2)
    ax_inset.plot(s_grid[peak_mask], pdf_exact[peak_mask], 'k-', linewidth=2)
    ax_inset.set_title('Peak region (0.6-1.4)', fontsize=10)
    ax_inset.grid(True, alpha=0.3)

    # --- Plot 2: Residuals ---
    residual_exact = hist_counts - pdf_exact_at_bins
    residual_wigner = hist_counts - pdf_wigner_at_bins

    ax2.plot(bin_centers, residual_wigner, 'r--', linewidth=1.5, alpha=0.7,
             label=f'Empirical - Wigner (MAD={mad_wigner:.4f})')
    ax2.plot(bin_centers, residual_exact, 'k-', linewidth=1.5,
             label=f'Empirical - Exact GUE (MAD={mad_exact:.4f})')
    ax2.axhline(0, color='gray', linewidth=0.5)
    ax2.set_xlabel('Normalized spacing s', fontsize=13)
    ax2.set_ylabel('Residual', fontsize=13)
    ax2.set_title('Residuals: Empirical minus Theoretical', fontsize=14)
    ax2.legend(fontsize=11)
    ax2.set_xlim(0, s_max)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = os.path.join(PLOTS_DIR, "exact_gue_comparison.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to {outpath}")

    # --- Summary ---
    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Exact GUE variance:                    {var_exact:.6f}")
    print(f"  Exact GUE mean:                        {mean_exact:.6f}")
    print(f"  Empirical spacing variance:            {np.var(spacings):.6f}")
    print(f"  Empirical spacing mean:                {np.mean(spacings):.6f}")
    print(f"  MAD from Wigner surmise:               {mad_wigner:.6f}")
    print(f"  MAD from exact GUE:                    {mad_exact:.6f}")
    print(f"  Exact GUE is {mad_wigner/mad_exact:.1f}x better fit than Wigner")


if __name__ == "__main__":
    main()

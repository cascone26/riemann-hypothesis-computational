"""
Wu-Sprung Potential Reconstruction

Wu & Sprung (1993) showed that a 1D Schrodinger operator
  H = -d^2/dx^2 + V(x) on [0, infinity) with psi(0) = 0
can be constructed from the Riemann zeros using WKB inversion.

The potential V_WS(x) is the unique 1D potential whose eigenvalues
match the zeta zeros. It has fractal structure.

This script:
1. Reconstructs V_WS(x) from the zero data
2. Verifies the eigenvalues of H match the input zeros
3. Analyzes the fractal dimension of V_WS
"""

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.linalg import eigvalsh
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)


def load_zeros(count=1000):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def N_smooth(E):
    """Smooth zero counting function."""
    if E < 2:
        return 0.0
    return (E / (2*np.pi)) * np.log(E / (2*np.pi)) - E / (2*np.pi) + 7/8


def reconstruct_potential_wkb(zeros, n_points=2000):
    """
    WKB inversion: reconstruct V(x) from the eigenvalues E_n = t_n^2
    (we work with E = t^2 since H = -d^2/dx^2 + V has eigenvalues E_n,
    and the zeta zeros t_n relate to eigenvalues).

    Actually, in the Wu-Sprung setup, the eigenvalues ARE the t_n directly
    (not t_n^2). The Hamiltonian is H = -d^2/dx^2 + V(x) with eigenvalues
    matching the imaginary parts of the zeros.

    The WKB quantization condition:
      integral_0^{x_turn(E)} sqrt(E - V(x)) dx = (n - 1/4) * pi

    for the n-th eigenvalue E_n, where x_turn(E) is the classical turning point.

    Inverting: x(V) = (1/pi) * dN/dE * integral... (Gel'fand-Levitan)

    Simpler approach: use the smooth counting function to get x(V):
      x(V) = integral_0^V (dN/dE) / sqrt(V - E) dE / pi

    This is Abel's integral equation inversion.
    """
    V_min = zeros[0] * 0.9  # below first zero
    V_max = zeros[-1] * 1.1  # above last zero

    V_values = np.linspace(V_min, V_max, n_points)
    x_values = np.zeros(n_points)

    for i, V in enumerate(V_values):
        if V <= zeros[0] * 0.5:
            x_values[i] = 0
            continue

        # x(V) = (1/pi) * integral_V0^V N'(E) / sqrt(V - E) dE
        # where N'(E) = dN/dE is the density of states
        def integrand(E):
            if E < 2 or V - E < 1e-10:
                return 0
            # N'(E) = (1/2pi) * log(E/2pi)
            density = np.log(E / (2*np.pi)) / (2*np.pi)
            return density / np.sqrt(V - E)

        try:
            result, _ = quad(integrand, max(2, V_min), V - 0.001,
                            limit=100, epsabs=1e-8)
            x_values[i] = result / np.pi
        except:
            x_values[i] = x_values[i-1] if i > 0 else 0

    return V_values, x_values


def verify_eigenvalues(V_values, x_values, zeros, n_verify=30):
    """
    Given V(x), solve the Schrodinger equation numerically
    and compare eigenvalues to the original zeros.
    """
    # Interpolate V(x) onto a regular grid
    x_max = x_values[-1]
    N = 2000
    dx = x_max / N
    x_grid = np.linspace(dx, x_max, N)

    # Interpolate V onto regular grid
    V_interp = np.interp(x_grid, x_values, V_values)

    # Build Hamiltonian matrix: H = -d^2/dx^2 + V(x)
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = 2.0 / dx**2 + V_interp[i]
        if i > 0:
            H[i, i-1] = -1.0 / dx**2
        if i < N-1:
            H[i, i+1] = -1.0 / dx**2

    # Compute eigenvalues
    eigenvalues = eigvalsh(H)
    eigenvalues = eigenvalues[:n_verify]

    # Compare
    n = min(len(eigenvalues), len(zeros), n_verify)
    target = zeros[:n]

    print(f"\n  Eigenvalue verification (N_grid={N}):")
    print(f"  {'#':>3} {'Computed':>12} {'Target Zero':>12} {'Rel Error':>10}")
    print(f"  {'-'*40}")

    for i in range(min(n, 15)):
        rel_err = abs(eigenvalues[i] - target[i]) / target[i]
        print(f"  {i+1:>3} {eigenvalues[i]:>12.4f} {target[i]:>12.4f} {rel_err:>10.6f}")

    mean_rel = np.mean(np.abs(eigenvalues[:n] - target[:n]) / target[:n])
    print(f"\n  Mean relative error: {mean_rel:.6f}")

    return eigenvalues[:n], mean_rel


def analyze_fractal_dimension(V_values, x_values):
    """
    Wu-Sprung found the potential has fractal dimension d ≈ 1.5.
    Estimate the box-counting dimension of V(x).
    """
    # Use the derivative dV/dx as a signal
    # Remove NaN/Inf
    valid = (x_values > 0) & np.isfinite(V_values)
    x = x_values[valid]
    V = V_values[valid]

    if len(x) < 100:
        print("  Not enough points for fractal analysis")
        return None

    # Box counting: how many boxes of size epsilon does the curve visit?
    def box_count(x, y, epsilon):
        x_boxes = np.floor(x / epsilon).astype(int)
        y_boxes = np.floor(y / epsilon).astype(int)
        boxes = set(zip(x_boxes, y_boxes))
        return len(boxes)

    # Normalize to [0,1] x [0,1]
    x_norm = (x - x.min()) / (x.max() - x.min() + 1e-10)
    V_norm = (V - V.min()) / (V.max() - V.min() + 1e-10)

    epsilons = np.logspace(-3, -0.5, 20)
    counts = [box_count(x_norm, V_norm, eps) for eps in epsilons]

    # Fit dimension
    log_eps = np.log(1/np.array(epsilons))
    log_counts = np.log(np.array(counts, dtype=float))

    valid_fit = np.isfinite(log_counts) & (np.array(counts) > 0)
    if sum(valid_fit) > 2:
        coeffs = np.polyfit(log_eps[valid_fit], log_counts[valid_fit], 1)
        fractal_dim = coeffs[0]
    else:
        fractal_dim = float('nan')

    print(f"\n  Fractal dimension analysis:")
    print(f"    Box-counting dimension: {fractal_dim:.4f}")
    print(f"    (Wu-Sprung found d ≈ 1.5)")

    return fractal_dim


def main():
    print("=" * 60)
    print("WU-SPRUNG POTENTIAL RECONSTRUCTION")
    print("=" * 60)

    zeros = load_zeros(count=500)
    print(f"Using {len(zeros)} zeros")
    print(f"Range: {zeros[0]:.4f} to {zeros[-1]:.4f}")

    # Reconstruct potential
    print("\nReconstructing V(x) via WKB inversion...")
    V_values, x_values = reconstruct_potential_wkb(zeros, n_points=3000)

    # Verify
    eigenvalues, mean_err = verify_eigenvalues(V_values, x_values, zeros, n_verify=30)

    # Fractal analysis
    fractal_dim = analyze_fractal_dimension(V_values, x_values)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # V(x)
    ax = axes[0, 0]
    valid = x_values > 0
    ax.plot(x_values[valid], V_values[valid], 'b-', linewidth=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('V(x)')
    ax.set_title('Wu-Sprung Potential')
    ax.grid(True, alpha=0.3)

    # V(x) zoomed — look for fractal structure
    ax = axes[0, 1]
    mid = len(x_values) // 2
    window = len(x_values) // 10
    ax.plot(x_values[mid-window:mid+window], V_values[mid-window:mid+window],
            'b-', linewidth=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('V(x)')
    ax.set_title('Wu-Sprung Potential (zoomed)')
    ax.grid(True, alpha=0.3)

    # Eigenvalue comparison
    ax = axes[1, 0]
    n = len(eigenvalues)
    ax.plot(range(1, n+1), zeros[:n], 'ro', markersize=4, label='Zeta zeros')
    ax.plot(range(1, n+1), eigenvalues, 'bx', markersize=6, label='H eigenvalues')
    ax.set_xlabel('n')
    ax.set_ylabel('Value')
    ax.set_title('Eigenvalues vs Zeta Zeros')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Relative errors
    ax = axes[1, 1]
    rel_errs = np.abs(eigenvalues - zeros[:n]) / zeros[:n]
    ax.plot(range(1, n+1), rel_errs, 'g.-')
    ax.set_xlabel('n')
    ax.set_ylabel('Relative error')
    ax.set_title('Eigenvalue Relative Errors')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    plt.suptitle('Wu-Sprung Potential Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'wu_sprung.png'), dpi=150)
    plt.close()
    print(f"\nPlots saved.")

    # Save results
    results = {
        "n_zeros_used": len(zeros),
        "mean_eigenvalue_error": float(mean_err),
        "fractal_dimension": float(fractal_dim) if fractal_dim else None,
        "n_potential_points": len(V_values),
    }
    with open(os.path.join(RESULTS_DIR, "wu_sprung.json"), "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()

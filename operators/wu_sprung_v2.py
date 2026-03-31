"""
Wu-Sprung Potential Reconstruction v2

Wu & Sprung (1993): the 1D Schrodinger equation
  -psi''(x) + V(x)*psi(x) = E_n * psi(x),  psi(0) = 0
has a unique potential V(x) whose eigenvalues E_n equal the Riemann zeros t_n.

The WKB inversion formula (Chadan-Sabatier):
  x(E) = (1/pi) * integral_{E_0}^{E} N(E') / sqrt(E - E') dE'

where N(E) is the zero counting function. This gives x as a function of E,
and we invert to get V(x) = E(x), the potential at the turning point.

The key insight: at the classical turning point, V(x_turn) = E, so the
curve x(E) IS the potential V(x) with axes swapped.
"""

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigvalsh
from scipy.interpolate import interp1d
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_zeros(count=500):
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
    """Smooth zero counting function N(E) ~ (E/2pi)*log(E/2pi) - E/2pi + 7/8"""
    if E <= 0:
        return 0.0
    return (E / (2*np.pi)) * np.log(E / (2*np.pi)) - E / (2*np.pi) + 7/8


def dN_dE(E):
    """Derivative of smooth counting function: density of states."""
    if E <= 0:
        return 0.0
    return np.log(E / (2*np.pi)) / (2*np.pi)


def reconstruct_x_of_E(n_points=2000):
    """
    WKB inversion: compute x(E) from the smooth counting function.

    x(E) = (1/pi) * integral_0^E [N'(E') / sqrt(E - E')] dE'

    This integral is the Abel transform of N'(E).
    x(E) gives the turning point position for energy E.
    Since V(x_turn) = E at the turning point, the curve (x(E), E)
    IS the potential V(x).
    """
    E_min = 5.0   # below first zero at ~14.1
    E_max = 500.0  # cover many zeros
    E_values = np.linspace(E_min, E_max, n_points)
    x_values = np.zeros(n_points)

    for i, E in enumerate(E_values):
        def integrand(Ep):
            if Ep <= 1 or E - Ep < 1e-10:
                return 0.0
            return dN_dE(Ep) / np.sqrt(E - Ep)

        try:
            # Split at singularity E' = E
            result, _ = quad(integrand, 1, E - 1e-6,
                           limit=200, epsabs=1e-10)
            x_values[i] = result / np.pi
        except:
            x_values[i] = x_values[i-1] if i > 0 else 0

    return E_values, x_values


def solve_schrodinger(V_func, x_max, N_grid=4000, n_eigs=50):
    """
    Solve -psi''(x) + V(x)*psi(x) = E*psi(x) on [0, x_max]
    with Dirichlet BCs: psi(0) = psi(x_max) = 0.
    """
    dx = x_max / N_grid
    x = np.linspace(dx, x_max - dx, N_grid - 1)  # interior points

    # Tridiagonal Hamiltonian
    diag = 2.0/dx**2 + V_func(x)
    off_diag = -np.ones(len(x) - 1) / dx**2

    H = np.diag(diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)

    eigenvalues = eigvalsh(H)
    return eigenvalues[:n_eigs]


def main():
    print("=" * 60)
    print("WU-SPRUNG POTENTIAL v2")
    print("=" * 60)

    zeros = load_zeros(200)
    print(f"Target: first {len(zeros)} zeros ({zeros[0]:.2f} to {zeros[-1]:.2f})")

    # Step 1: Reconstruct x(E)
    print("\nStep 1: Computing x(E) via WKB Abel inversion...")
    E_values, x_values = reconstruct_x_of_E(n_points=3000)

    # x(E) should be monotonically increasing
    # Clean up any non-monotonic parts
    for i in range(1, len(x_values)):
        if x_values[i] <= x_values[i-1]:
            x_values[i] = x_values[i-1] + 1e-6

    print(f"  x range: {x_values[0]:.4f} to {x_values[-1]:.4f}")
    print(f"  E range: {E_values[0]:.4f} to {E_values[-1]:.4f}")

    # Step 2: V(x) = E(x) = inverse of x(E)
    # Interpolate: given x, find E such that x(E) = x
    V_of_x = interp1d(x_values, E_values, kind='cubic',
                       bounds_error=False, fill_value='extrapolate')

    x_max = x_values[-1]
    print(f"  Turning point at E={E_values[-1]:.1f}: x = {x_max:.4f}")

    # Step 3: Solve the Schrodinger equation with this potential
    print("\nStep 2: Solving Schrodinger equation...")
    eigenvalues = solve_schrodinger(V_of_x, x_max, N_grid=5000, n_eigs=50)

    # Step 4: Compare
    n = min(len(eigenvalues), len(zeros), 30)
    target = zeros[:n]
    eigs = eigenvalues[:n]

    print(f"\nStep 3: Comparing eigenvalues to zeta zeros")
    print(f"  {'#':>3} {'Eigenvalue':>12} {'Zero':>12} {'Rel Error':>10}")
    print(f"  {'-'*40}")

    for i in range(min(n, 20)):
        rel_err = abs(eigs[i] - target[i]) / target[i]
        marker = "  " if rel_err < 0.05 else " *"
        print(f"  {i+1:>3} {eigs[i]:>12.4f} {target[i]:>12.4f} {rel_err:>10.6f}{marker}")

    rel_errors = np.abs(eigs[:n] - target[:n]) / target[:n]
    mean_err = np.mean(rel_errors)
    print(f"\n  Mean relative error (first {n}): {mean_err:.6f}")
    print(f"  Max relative error: {np.max(rel_errors):.6f}")

    # Step 5: Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # V(x) — the Wu-Sprung potential
    ax = axes[0, 0]
    x_plot = np.linspace(x_values[0], x_values[-1], 1000)
    V_plot = V_of_x(x_plot)
    ax.plot(x_plot, V_plot, 'b-', linewidth=1)
    # Mark the zeta zeros as horizontal lines at the eigenvalue heights
    for i, z in enumerate(zeros[:20]):
        ax.axhline(y=z, color='red', alpha=0.2, linewidth=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('V(x)')
    ax.set_title('Wu-Sprung Potential (WKB smooth part)')
    ax.grid(True, alpha=0.3)

    # x(E) — the turning point function
    ax = axes[0, 1]
    ax.plot(E_values, x_values, 'b-', linewidth=1)
    for z in zeros[:20]:
        ax.axvline(x=z, color='red', alpha=0.2, linewidth=0.5)
    ax.set_xlabel('E (energy / zero height)')
    ax.set_ylabel('x(E) (turning point)')
    ax.set_title('Classical Turning Points')
    ax.grid(True, alpha=0.3)

    # Eigenvalue comparison
    ax = axes[1, 0]
    ax.plot(range(1, n+1), target, 'ro', markersize=5, label='Zeta zeros', zorder=5)
    ax.plot(range(1, n+1), eigs[:n], 'bx', markersize=7, label='Schrodinger eigenvalues')
    ax.set_xlabel('n')
    ax.set_ylabel('Value')
    ax.set_title('Eigenvalues vs Zeta Zeros')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Relative errors
    ax = axes[1, 1]
    ax.plot(range(1, n+1), rel_errors, 'g.-')
    ax.set_xlabel('n')
    ax.set_ylabel('Relative error |E_n - t_n| / t_n')
    ax.set_title(f'Eigenvalue Errors (mean: {mean_err:.4f})')
    ax.grid(True, alpha=0.3)

    plt.suptitle('Wu-Sprung Potential — Smooth WKB Reconstruction', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'wu_sprung_v2.png'), dpi=150)
    plt.close()
    print(f"\nPlots saved to {PLOTS_DIR}/wu_sprung_v2.png")

    # Save
    results = {
        "n_zeros_target": int(len(zeros)),
        "n_eigenvalues": int(n),
        "mean_rel_error": float(mean_err),
        "max_rel_error": float(np.max(rel_errors)),
        "first_20": [
            {"n": i+1, "eigenvalue": float(eigs[i]), "zero": float(target[i]),
             "rel_error": float(rel_errors[i])}
            for i in range(min(20, n))
        ],
    }
    with open(os.path.join(RESULTS_DIR, "wu_sprung_v2.json"), "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()

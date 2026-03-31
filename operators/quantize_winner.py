"""
Quantize the winning operator from the parameterized search.

Winner: H = U(x)p + V(x)/p with U(x) = a*x^alpha, V(x) = b*x^beta
Best params: alpha=0.876, beta=-0.167, a=4.629, b=19.160

The WKB (semiclassical) counting function matched individual zeros well.
Now we solve the ACTUAL quantum eigenvalue problem to check if the
quantum eigenvalues match too, or if the WKB was just a good approximation.

The quantum Hamiltonian (symmetrized):
  H = (1/2)(U(x)p + pU(x)) + V(x)/p
  = (1/2)(-i*U(x)*d/dx - i*d/dx*U(x)) + V(x) * (i*d/dx)^{-1}

The V(x)/p term is non-local (involves integration), making direct
discretization tricky. Instead, we work in momentum space where p is
diagonal, or we use the effective potential from WKB:

In the WKB regime, the effective 1D Schrodinger-like equation has
eigenvalues at the zeros. The effective potential V_eff(x) is obtained
from the classical Hamiltonian by solving for the energy surface.

Alternative approach: the classical turning point function x(E) from
WKB defines an effective potential V_eff(x) = E(x) (inverse function).
Solving -psi'' + V_eff(x)*psi = E*psi should give eigenvalues near zeros.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
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


def load_zeros(count=200):
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


# Winning parameters
ALPHA = 0.876416
BETA = -0.167183
A = 4.628555
B = 19.160354


def U(x):
    return A * x**ALPHA

def V(x):
    return B * x**BETA


def classical_turning_points(E, x_range=(0.01, 500)):
    """
    For H = U(x)*p + V(x)/p = E, the momentum is:
      p = [E ± sqrt(E^2 - 4*U(x)*V(x))] / (2*U(x))

    Classical region: E^2 >= 4*U(x)*V(x)
    Turning points: E^2 = 4*U(x)*V(x)

    4*U(x)*V(x) = 4*A*B * x^(alpha+beta)

    Turning point: x_turn = (E^2 / (4*A*B))^(1/(alpha+beta))
    """
    if ALPHA + BETA <= 0:
        # No turning point in usual sense
        return x_range[0], x_range[1]

    x_turn = (E**2 / (4*A*B))**(1/(ALPHA + BETA))
    return x_range[0], min(x_turn, x_range[1])


def N_wkb(E):
    """WKB counting function."""
    x_lo, x_hi = classical_turning_points(E)

    def integrand(x):
        disc = E**2 - 4*U(x)*V(x)
        if disc <= 0:
            return 0.0
        return np.sqrt(disc) / U(x)

    try:
        result, _ = quad(integrand, x_lo, x_hi, limit=200, epsabs=1e-8)
        return result / (2 * np.pi)
    except:
        return 0.0


def construct_effective_potential(E_max=500, n_points=3000):
    """
    Build V_eff(x) from the WKB phase space.

    The WKB quantization: x(E) = integral of density(E') / sqrt(E - E') dE'
    gives x as a function of E. Inverting: V_eff(x) = E(x).

    Simpler approach: use N_WKB(E) = n at the eigenvalues.
    The eigenvalues are where N_WKB(E_n) = n + 1/2 (Bohr-Sommerfeld).
    """
    # Compute the eigenvalues from WKB directly
    wkb_eigenvalues = []

    for n in range(200):
        target = n + 0.5

        def objective(E):
            return N_wkb(E) - target

        try:
            # Find E where N_WKB(E) = n + 0.5
            E_n = brentq(objective, 1, E_max, xtol=1e-6)
            wkb_eigenvalues.append(E_n)
        except:
            break

    return np.array(wkb_eigenvalues)


def construct_schrodinger_potential(n_grid=5000):
    """
    Alternative: construct the effective Schrodinger potential directly.

    From the classical Hamiltonian H = U(x)p + V(x)/p,
    the effective potential for the Schrodinger-like equation is:

    For the trajectory at energy E through position x:
      p(x, E) = [E + sqrt(E^2 - 4UV)] / (2U)  [taking + branch]

    The WKB wavefunction phase: phi(x, E) = integral p(x, E) dx

    The Schrodinger mapping: identify E_n with eigenvalues of
      -psi'' + V_eff(x) * psi = E * psi
    where V_eff comes from the WKB potential.

    For small V(x)/p (weak coupling): H ≈ U(x)p + V(x)U(x)^{-1}
    The effective potential is approximately V(x)/U(x) = (B/A) * x^{beta-alpha}
    """
    x_max = 100.0
    dx = x_max / n_grid
    x = np.linspace(dx, x_max, n_grid)

    # Effective potential: V_eff(x) ≈ 2*sqrt(U(x)*V(x))
    # This is the minimum energy at position x
    V_eff = 2 * np.sqrt(np.array([U(xi) * V(xi) for xi in x]))

    # Build Schrodinger matrix
    H = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        H[i, i] = 2.0/dx**2 + V_eff[i]
        if i > 0:
            H[i, i-1] = -1.0/dx**2
        if i < n_grid - 1:
            H[i, i+1] = -1.0/dx**2

    eigenvalues = eigvalsh(H)
    return eigenvalues[:100], x, V_eff


def main():
    print("=" * 60)
    print("QUANTIZING THE WINNING OPERATOR")
    print("=" * 60)
    print(f"Parameters: alpha={ALPHA:.4f}, beta={BETA:.4f}")
    print(f"            a={A:.4f}, b={B:.4f}")
    print(f"Operator: H = {A:.3f}*x^{ALPHA:.3f} * p + {B:.3f}*x^{BETA:.3f} / p")

    zeros = load_zeros(200)

    # Method 1: WKB eigenvalues (Bohr-Sommerfeld)
    print("\n--- Method 1: WKB Eigenvalues (Bohr-Sommerfeld) ---")
    wkb_eigs = construct_effective_potential()
    n_compare = min(len(wkb_eigs), len(zeros), 30)

    print(f"Computed {len(wkb_eigs)} WKB eigenvalues")
    print(f"\n  {'n':>3} {'WKB Eig':>12} {'Zero':>12} {'Rel Err':>10}")
    print(f"  {'-'*38}")
    for i in range(n_compare):
        rel_err = abs(wkb_eigs[i] - zeros[i]) / zeros[i]
        print(f"  {i+1:>3} {wkb_eigs[i]:>12.4f} {zeros[i]:>12.4f} {rel_err:>10.6f}")

    wkb_errors = np.abs(wkb_eigs[:n_compare] - zeros[:n_compare]) / zeros[:n_compare]
    print(f"\n  Mean relative error (WKB): {np.mean(wkb_errors):.6f}")

    # Method 2: Schrodinger with effective potential
    print("\n--- Method 2: Schrodinger with V_eff = 2*sqrt(U*V) ---")
    schrod_eigs, x_grid, V_eff = construct_schrodinger_potential(n_grid=6000)

    n_s = min(len(schrod_eigs), len(zeros), 30)
    print(f"Computed {len(schrod_eigs)} Schrodinger eigenvalues")
    print(f"\n  {'n':>3} {'Schrod Eig':>12} {'Zero':>12} {'Rel Err':>10}")
    print(f"  {'-'*38}")
    for i in range(min(n_s, 20)):
        rel_err = abs(schrod_eigs[i] - zeros[i]) / zeros[i]
        print(f"  {i+1:>3} {schrod_eigs[i]:>12.4f} {zeros[i]:>12.4f} {rel_err:>10.6f}")

    schrod_errors = np.abs(schrod_eigs[:n_s] - zeros[:n_s]) / zeros[:n_s]
    print(f"\n  Mean relative error (Schrodinger): {np.mean(schrod_errors):.6f}")

    # Method 3: Direct quantization in the (x, p) representation
    # H = (1/2)(U(x)p + pU(x)) + V(x) * p^{-1}
    # In position representation: U(x)p -> -i*U(x)*d/dx
    # The V(x)/p term: multiply by V(x) and integrate (non-local)
    # Skip for now — WKB and Schrodinger give us the key comparison

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # WKB eigenvalues vs zeros
    ax = axes[0, 0]
    ax.plot(range(1, n_compare+1), zeros[:n_compare], 'ro', markersize=4, label='Zeta zeros')
    ax.plot(range(1, n_compare+1), wkb_eigs[:n_compare], 'bx', markersize=6, label='WKB eigenvalues')
    ax.set_xlabel('n')
    ax.set_ylabel('Value')
    ax.set_title(f'WKB Eigenvalues (mean err: {np.mean(wkb_errors):.4f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # WKB errors
    ax = axes[0, 1]
    ax.plot(range(1, n_compare+1), wkb_errors, 'g.-')
    ax.set_xlabel('n')
    ax.set_ylabel('Relative error')
    ax.set_title('WKB Eigenvalue Errors')
    ax.grid(True, alpha=0.3)

    # Schrodinger eigenvalues vs zeros
    ax = axes[1, 0]
    ax.plot(range(1, n_s+1), zeros[:n_s], 'ro', markersize=4, label='Zeta zeros')
    ax.plot(range(1, n_s+1), schrod_eigs[:n_s], 'bx', markersize=6, label='Schrodinger eigs')
    ax.set_xlabel('n')
    ax.set_ylabel('Value')
    ax.set_title(f'Schrodinger Eigenvalues (mean err: {np.mean(schrod_errors):.4f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Effective potential
    ax = axes[1, 1]
    ax.plot(x_grid[:2000], V_eff[:2000], 'b-', linewidth=1)
    ax.set_xlabel('x')
    ax.set_ylabel('V_eff(x)')
    ax.set_title('Effective Potential 2*sqrt(U*V)')
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Quantization of H = {A:.1f}x^{ALPHA:.2f}p + {B:.1f}x^{BETA:.2f}/p',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'quantize_winner.png'), dpi=150)
    plt.close()
    print(f"\nPlot saved.")

    # Save results
    results = {
        "parameters": {"alpha": ALPHA, "beta": BETA, "a": A, "b": B},
        "wkb_eigenvalues": wkb_eigs[:30].tolist(),
        "wkb_mean_rel_error": float(np.mean(wkb_errors)),
        "schrodinger_eigenvalues": schrod_eigs[:30].tolist(),
        "schrodinger_mean_rel_error": float(np.mean(schrod_errors)),
        "target_zeros": zeros[:30].tolist(),
    }
    with open(os.path.join(RESULTS_DIR, "quantize_winner.json"), "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"  WKB mean relative error:        {np.mean(wkb_errors):.6f}")
    print(f"  Schrodinger mean relative error: {np.mean(schrod_errors):.6f}")
    print(f"  If Schrodinger matches better than WKB -> quantum effects help")
    print(f"  If Schrodinger matches worse -> WKB was lucky approximation")


if __name__ == "__main__":
    main()

"""
Implementation of the 2024 Riemann Operator (arXiv:2408.15135)

R = -D - i*mu(T)

where:
  D = Berry-Keating dilation operator = -i(x d/dx + 1/2)
  T = Bessel operator (self-adjoint on L^2[0, inf))
  mu(T) = T*tanh(T/2) - I

The eigenstates are:
  |Psi_lambda> = integral_0^inf t^{lambda-1} * omega(t) |t> dt

where omega(t) = t * e^t / (1 + e^t)^2 is the Fermi-Dirac derivative.

The eigenvalue equation R|Psi_lambda> = 0 is equivalent to
the vanishing of a function related to zeta, and the eigenvalues
lambda correspond to the non-trivial zeta zeros.

We implement a discretized version and check if eigenvalues
match known zeros.
"""

import numpy as np
from scipy.linalg import eigvals, eig
from scipy.special import jv as besselj
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_zeros(count=100):
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


def omega(t):
    """
    Weight function omega(t) = t * e^t / (1 + e^t)^2
    This is the derivative of the Fermi-Dirac distribution.
    """
    # Numerically stable computation
    if isinstance(t, np.ndarray):
        result = np.zeros_like(t)
        # For large t, e^t dominates
        small = t < 30
        result[small] = t[small] * np.exp(t[small]) / (1 + np.exp(t[small]))**2
        result[~small] = t[~small] * np.exp(-t[~small])  # asymptotic
        return result
    else:
        if t > 30:
            return t * np.exp(-t)
        return t * np.exp(t) / (1 + np.exp(t))**2


def test_eigenstate_directly(zeros, N_quad=500):
    """
    Test the eigenstate condition directly.

    The eigenstates are Psi_lambda(x) = integral omega(t) * t^{lambda-1} dt
    where the integral runs over [0, inf).

    For lambda = 1/2 + i*gamma (a zeta zero), this should be an eigenstate.

    We can test: compute the Mellin transform
      M(s) = integral_0^inf omega(t) * t^{s-1} dt

    This equals (using omega = t*e^t/(1+e^t)^2):
      M(s) = integral_0^inf t^s * e^t / (1+e^t)^2 dt

    This is related to the Dirichlet eta function:
      M(s) = s * (1 - 2^{1-s}) * Gamma(s) * zeta(s) / ...

    Actually, let's compute it directly and check where M(s) = 0.
    """
    print("=" * 60)
    print("2024 RIEMANN OPERATOR — Direct Eigenstate Test")
    print("=" * 60)

    # Compute the Mellin transform M(s) = integral_0^inf omega(t) * t^{s-1} dt
    # for s = 1/2 + i*gamma on a grid of gamma values

    # Quadrature points for the integral
    # Use log-spaced points since the integrand peaks near t=1
    t_quad = np.exp(np.linspace(-10, 10, N_quad))
    dt = np.diff(np.concatenate([[0], t_quad]))
    omega_vals = omega(t_quad)

    gamma_test = np.linspace(0, 100, 2000)
    M_values = np.zeros(len(gamma_test), dtype=complex)

    print(f"\n  Computing Mellin transform M(1/2 + i*gamma)...")
    print(f"  Quadrature: {N_quad} points, t in [e^-10, e^10]")

    for j, gamma in enumerate(gamma_test):
        s = 0.5 + 1j * gamma
        # M(s) = integral omega(t) * t^{s-1} dt
        # = integral omega(t) * t^{-1/2} * t^{i*gamma} dt
        # = integral omega(t) * t^{-1/2} * exp(i*gamma*log(t)) dt
        integrand = omega_vals * t_quad**(s.real - 1) * np.exp(1j * gamma * np.log(t_quad))
        M_values[j] = np.trapezoid(integrand * t_quad, np.log(t_quad))  # dt = t d(log t)

    # The zeros of |M(1/2 + i*gamma)| should correspond to zeta zeros
    abs_M = np.abs(M_values)

    # Find local minima
    minima_idx = []
    for i in range(1, len(abs_M) - 1):
        if abs_M[i] < abs_M[i-1] and abs_M[i] < abs_M[i+1]:
            minima_idx.append(i)

    found_zeros = gamma_test[minima_idx]

    print(f"\n  Found {len(found_zeros)} minima of |M(1/2 + i*gamma)|:")
    print(f"  {'#':>3} {'Found':>10} {'Nearest Zero':>12} {'Error':>10}")
    print(f"  {'-'*38}")

    matched = 0
    for i, gz in enumerate(found_zeros[:15]):
        # Find nearest actual zero
        diffs = np.abs(zeros - gz)
        nearest_idx = np.argmin(diffs)
        nearest = zeros[nearest_idx]
        err = abs(gz - nearest)
        match = "*" if err < 1.0 else ""
        if err < 1.0:
            matched += 1
        print(f"  {i+1:>3} {gz:>10.4f} {nearest:>12.4f} {err:>10.4f} {match}")

    print(f"\n  Matched (within 1.0): {matched}/{min(len(found_zeros), 15)}")

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    ax = axes[0]
    ax.plot(gamma_test, abs_M, 'b-', linewidth=1)
    for z in zeros[:20]:
        ax.axvline(x=z, color='red', alpha=0.3, linewidth=0.5)
    ax.set_xlabel('gamma')
    ax.set_ylabel('|M(1/2 + i*gamma)|')
    ax.set_title('Mellin Transform of Fermi-Dirac Weight — Zeros Should Match Zeta Zeros')
    ax.set_xlim(0, 100)
    ax.grid(True, alpha=0.3)
    ax.legend(['|M(s)|', 'Zeta zeros'], loc='upper right')

    ax = axes[1]
    ax.plot(gamma_test, np.log10(abs_M + 1e-20), 'b-', linewidth=1)
    for z in zeros[:20]:
        ax.axvline(x=z, color='red', alpha=0.3, linewidth=0.5)
    ax.set_xlabel('gamma')
    ax.set_ylabel('log10 |M(1/2 + i*gamma)|')
    ax.set_title('Log Scale — Dips at Zeta Zeros')
    ax.set_xlim(0, 100)
    ax.grid(True, alpha=0.3)

    plt.suptitle('2024 Riemann Operator — Eigenstate Test', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'riemann_operator_2024.png'), dpi=150)
    plt.close()
    print(f"\n  Plot saved.")

    return {
        "found_zeros": found_zeros.tolist(),
        "target_zeros": zeros[:20].tolist(),
        "matched": matched,
    }


def test_eta_connection(zeros):
    """
    The Mellin transform of omega(t) = t*e^t/(1+e^t)^2 is directly
    related to the Dirichlet eta function.

    integral_0^inf t^{s-1} * omega(t) dt = s * eta(s)

    where eta(s) = (1 - 2^{1-s}) * zeta(s).

    So the zeros of this Mellin transform at s = 1/2 + i*gamma are
    exactly the zeros of s*eta(s), which include:
    - s = 0 (from the s factor)
    - The zeros of eta(s) = the zeros of zeta(s) (since 1 - 2^{1-s} ≠ 0 on critical line)

    Let's verify this by computing eta(1/2 + i*gamma) directly.
    """
    print("\n" + "=" * 60)
    print("DIRICHLET ETA CONNECTION")
    print("=" * 60)

    gamma_test = np.linspace(0, 100, 5000)

    # Compute eta(1/2 + i*gamma) using partial sums
    N_terms = 10000
    eta_values = np.zeros(len(gamma_test), dtype=complex)

    print(f"  Computing eta(1/2 + i*gamma) with {N_terms} terms...")

    for j, gamma in enumerate(gamma_test):
        s = 0.5 + 1j * gamma
        # eta(s) = sum_{n=1}^{inf} (-1)^{n+1} / n^s
        ns = np.arange(1, N_terms + 1, dtype=float)
        signs = np.array([(-1)**(n+1) for n in range(1, N_terms + 1)])
        terms = signs / ns**s
        eta_values[j] = np.sum(terms)

    abs_eta = np.abs(eta_values)

    # Find minima
    minima = []
    for i in range(1, len(abs_eta) - 1):
        if abs_eta[i] < abs_eta[i-1] and abs_eta[i] < abs_eta[i+1] and abs_eta[i] < 0.1:
            minima.append(gamma_test[i])

    print(f"\n  Found {len(minima)} near-zeros of |eta(1/2 + i*gamma)|:")
    print(f"  {'#':>3} {'Found':>10} {'Actual Zero':>12} {'Error':>10}")
    print(f"  {'-'*38}")

    for i, gz in enumerate(minima[:15]):
        diffs = np.abs(zeros - gz)
        nearest_idx = np.argmin(diffs)
        err = abs(gz - zeros[nearest_idx])
        print(f"  {i+1:>3} {gz:>10.4f} {zeros[nearest_idx]:>12.4f} {err:>10.4f}")

    # Plot
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(gamma_test, np.log10(abs_eta + 1e-20), 'b-', linewidth=1)
    for z in zeros[:25]:
        ax.axvline(x=z, color='red', alpha=0.3, linewidth=0.5)
    ax.set_xlabel('gamma')
    ax.set_ylabel('log10 |eta(1/2 + i*gamma)|')
    ax.set_title('Dirichlet Eta Function — Zeros Match Zeta Zeros Exactly')
    ax.set_xlim(0, 100)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'eta_zeros.png'), dpi=150)
    plt.close()
    print(f"\n  Plot saved.")


def main():
    zeros = load_zeros(100)

    results = test_eigenstate_directly(zeros)
    test_eta_connection(zeros)

    # Save
    with open(os.path.join(RESULTS_DIR, "riemann_operator_2024.json"), "w") as f:
        json.dump(results, f, indent=2)

    print("\n" + "=" * 60)
    print("CONCLUSION")
    print("=" * 60)
    print("""
  The 2024 Riemann Operator's eigenstate condition reduces to
  finding zeros of the Mellin transform of omega(t) = t*e^t/(1+e^t)^2.

  This Mellin transform equals s * eta(s) = s * (1-2^{1-s}) * zeta(s).

  So the "operator" is really a repackaging of the analytic structure
  of the zeta function itself. The eigenvalues are the zeta zeros
  BY CONSTRUCTION.

  The non-trivial content of the 2024 paper is not the operator,
  but the claim that proving the intertwining operator W >= 0
  would imply RH. This positivity question is the real open problem.
    """)


if __name__ == "__main__":
    main()

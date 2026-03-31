"""
Nyman-Beurling Distance to RH

RH is true iff the characteristic function chi_{(0,1]} can be approximated
in L^2(0,1) by functions of the form:
  f(t) = sum_k c_k * {theta_k / t}

where {x} = x - floor(x) is the fractional part.

The L^2 distance d_N = inf ||chi - f||_2 over all such f with N terms
should converge to 0 as N -> infinity iff RH is true.

Baez-Duarte showed we can restrict theta_k = 1/k (integers only).
His distance sequence d_N is a direct numerical measure of
"how close RH is to being true."

This script computes d_N for increasing N.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad
import json
import os

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")


def fractional_part(x):
    """Fractional part {x} = x - floor(x)."""
    return x - np.floor(x)


def inner_product_rho_rho(theta_j, theta_k):
    """
    <rho(theta_j/t), rho(theta_k/t)>_{L^2(0,1)}

    Analytical formula (Baez-Duarte):
    = 1 - (1/2)(theta_j + theta_k) + theta_j * theta_k * (1 - log(max(theta_j,theta_k)))

    But this only works for theta_j, theta_k in (0,1].
    For general theta_j = 1/j, theta_k = 1/k:

    <rho(1/(jt)), rho(1/(kt))>_{L^2(0,1)}
    = integral_0^1 {1/(jt)} * {1/(kt)} dt

    We compute this numerically for correctness.
    """
    result, _ = quad(lambda t: fractional_part(theta_j/t) * fractional_part(theta_k/t),
                     1e-10, 1, limit=200)
    return result


def inner_product_rho_chi(theta):
    """
    <rho(theta/t), chi_{(0,1]}>_{L^2(0,1)}
    = integral_0^1 {theta/t} dt
    """
    result, _ = quad(lambda t: fractional_part(theta/t), 1e-10, 1, limit=200)
    return result


def compute_gram_matrix(thetas):
    """Build the Gram matrix G[j,k] = <rho_j, rho_k> and vector b[j] = <rho_j, chi>."""
    N = len(thetas)
    G = np.zeros((N, N))
    b = np.zeros(N)

    for j in range(N):
        b[j] = inner_product_rho_chi(thetas[j])
        for k in range(j, N):
            G[j, k] = inner_product_rho_rho(thetas[j], thetas[k])
            G[k, j] = G[j, k]

    return G, b


def nyman_beurling_distance(N_max=30):
    """
    Compute d_N^2 = 1 - b^T G^{-1} b

    where G is the Gram matrix and b is the projection vector.
    d_N^2 = ||chi||^2 - max correlation = 1 - b^T G^{-1} b

    Using Baez-Duarte's integer restriction: theta_k = 1/k for k = 1, 2, ..., N
    """
    print("=" * 60)
    print("NYMAN-BEURLING DISTANCE TO RH")
    print("=" * 60)
    print(f"Computing d_N for N = 1 to {N_max}")
    print(f"Using Baez-Duarte restriction: theta_k = 1/k")
    print()

    # ||chi||^2 = integral_0^1 1 dt = 1
    chi_norm_sq = 1.0

    results = []

    for N in range(1, N_max + 1):
        thetas = [1.0 / k for k in range(1, N + 1)]

        G, b = compute_gram_matrix(thetas)

        try:
            # Add small regularization for numerical stability
            G_reg = G + 1e-12 * np.eye(N)
            coeffs = np.linalg.solve(G_reg, b)
            d_sq = chi_norm_sq - np.dot(b, coeffs)
            d_sq = max(d_sq, 0)  # numerical floor at 0
            d = np.sqrt(d_sq)
        except np.linalg.LinAlgError:
            d = float('nan')
            d_sq = float('nan')
            coeffs = np.zeros(N)

        results.append({
            "N": N,
            "d_squared": float(d_sq),
            "d": float(d),
            "coefficients": coeffs.tolist()[:5],  # first 5 only
        })

        print(f"  N={N:>3d}: d_N^2 = {d_sq:.10f}, d_N = {d:.10f}")

    # Check convergence rate
    print("\n  Convergence analysis:")
    for i in range(1, len(results)):
        if results[i-1]["d"] > 0 and results[i]["d"] > 0:
            ratio = results[i]["d"] / results[i-1]["d"]
            print(f"    d_{results[i]['N']} / d_{results[i-1]['N']} = {ratio:.6f}")

    # Plot
    Ns = [r["N"] for r in results]
    ds = [r["d"] for r in results]
    d_sqs = [r["d_squared"] for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.plot(Ns, ds, 'bo-', markersize=4)
    ax1.set_xlabel('N (number of terms)')
    ax1.set_ylabel('d_N')
    ax1.set_title('Nyman-Beurling Distance to RH')
    ax1.grid(True, alpha=0.3)

    ax2.semilogy(Ns, [max(d, 1e-15) for d in d_sqs], 'ro-', markersize=4)
    ax2.set_xlabel('N (number of terms)')
    ax2.set_ylabel('d_N^2 (log scale)')
    ax2.set_title('Distance Squared (RH true iff -> 0)')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'nyman_beurling_distance.png'), dpi=150)
    plt.close()
    print(f"\n  Plot saved.")

    # Save
    with open(os.path.join(RESULTS_DIR, "nyman_beurling.json"), "w") as f:
        json.dump(results, f, indent=2)

    return results


if __name__ == "__main__":
    nyman_beurling_distance(N_max=25)

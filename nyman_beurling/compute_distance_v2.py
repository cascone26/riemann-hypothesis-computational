"""
Nyman-Beurling Distance — v2 (incremental, N up to 200)

Fix: build the Gram matrix ONCE up to N_max, then solve N×N subsystems.
Original code rebuilt the full matrix from scratch at each N — O(N³) total integrals.
This version is O(N²) total integrals, making N=200 feasible.

RH ↔ d_N → 0. We want to see the convergence rate clearly.
"""

import numpy as np
from scipy.integrate import quad
import json
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR   = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
N_MAX = 200


def frac(x):
    return x - np.floor(x)


def G_entry(j, k):
    """∫₀¹ {1/(j·t)} {1/(k·t)} dt  — numerically via quad."""
    result, _ = quad(
        lambda t: frac(1.0 / (j * t)) * frac(1.0 / (k * t)),
        1e-10, 1.0, limit=200, epsabs=1e-10, epsrel=1e-10
    )
    return result


def b_entry(k):
    """∫₀¹ {1/(k·t)} dt  — numerically via quad."""
    result, _ = quad(
        lambda t: frac(1.0 / (k * t)),
        1e-10, 1.0, limit=200, epsabs=1e-10, epsrel=1e-10
    )
    return result


def build_gram_matrix(N_max):
    """
    Build (N_max × N_max) Gram matrix G and length-N_max vector b.
    Uses symmetry: G[j,k] = G[k,j]. O(N_max²/2) integrals.
    """
    G = np.zeros((N_max, N_max))
    b = np.zeros(N_max)

    t0 = time.time()
    total = N_max * (N_max + 1) // 2 + N_max
    done = 0

    # b vector
    for k in range(N_max):
        b[k] = b_entry(k + 1)  # k+1 because we use 1-indexed k
        done += 1
        if done % 20 == 0:
            elapsed = time.time() - t0
            rate = done / elapsed
            remaining = (total - done) / rate
            print(f"  [{done}/{total}] b entries... ({remaining:.0f}s remaining)", flush=True)

    # G matrix (upper triangle)
    for j in range(N_max):
        for k in range(j, N_max):
            G[j, k] = G_entry(j + 1, k + 1)
            G[k, j] = G[j, k]
            done += 1
            if done % 200 == 0:
                elapsed = time.time() - t0
                rate = done / elapsed
                remaining = (total - done) / rate
                print(f"  [{done}/{total}] G entries (j={j+1},k={k+1})... ({remaining:.0f}s remaining)", flush=True)

    elapsed = time.time() - t0
    print(f"  Done. Total time: {elapsed:.1f}s")
    return G, b


def solve_distances(G, b, N_max):
    """
    For each N from 1 to N_max, solve the N×N system and compute d_N.
    """
    results = []
    for N in range(1, N_max + 1):
        Gn = G[:N, :N] + 1e-13 * np.eye(N)
        bn = b[:N]
        try:
            coeffs = np.linalg.solve(Gn, bn)
            d_sq = max(1.0 - np.dot(bn, coeffs), 0.0)
            d = np.sqrt(d_sq)
        except np.linalg.LinAlgError:
            d_sq = float('nan')
            d = float('nan')
        results.append({"N": N, "d_sq": float(d_sq), "d": float(d)})
        if N % 10 == 0 or N <= 5:
            print(f"  N={N:>3d}: d_N = {d:.8f}  (d_N² = {d_sq:.2e})", flush=True)
    return results


def fit_convergence(results):
    Ns  = np.array([r["N"] for r in results if r["d"] > 1e-14])
    ds  = np.array([r["d"] for r in results if r["d"] > 1e-14])
    log_N = np.log(Ns)
    log_d = np.log(ds)

    # Power law: d_N ~ A * N^(-alpha)
    coeffs = np.polyfit(log_N, log_d, 1)
    alpha = -coeffs[0]
    A = np.exp(coeffs[1])

    # Also fit last half only (asymptotic rate)
    half = len(Ns) // 2
    coeffs2 = np.polyfit(log_N[half:], log_d[half:], 1)
    alpha2 = -coeffs2[0]
    A2 = np.exp(coeffs2[1])

    print(f"\n  Power law fit (all N):     d_N ~ {A:.4f} * N^(-{alpha:.4f})")
    print(f"  Power law fit (N>{Ns[half]:.0f}):  d_N ~ {A2:.4f} * N^(-{alpha2:.4f})")
    return alpha, A, alpha2, A2


def main():
    print("NYMAN-BEURLING DISTANCE v2 — N up to 200")
    print("=" * 60)
    print(f"Building {N_MAX}×{N_MAX} Gram matrix...")
    print()

    G, b = build_gram_matrix(N_MAX)

    print("\nSolving N×N subsystems for d_N...")
    results = solve_distances(G, b, N_MAX)

    print("\nConvergence analysis:")
    alpha, A, alpha2, A2 = fit_convergence(results)

    # Print summary at key N values
    print("\nSummary at key N:")
    key_Ns = [1, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200]
    for r in results:
        if r["N"] in key_Ns:
            print(f"  N={r['N']:>3d}: d_N = {r['d']:.8f}  d_N² = {r['d_sq']:.2e}")

    # Save
    with open(os.path.join(RESULTS_DIR, "nyman_beurling_200.json"), "w") as f:
        json.dump({
            "N_max": N_MAX,
            "power_law_alpha_all": float(alpha),
            "power_law_A_all": float(A),
            "power_law_alpha_asymptotic": float(alpha2),
            "power_law_A_asymptotic": float(A2),
            "results": results
        }, f, indent=2)

    # Plot
    Ns_plot = [r["N"] for r in results]
    ds_plot = [r["d"] for r in results]
    dsq_plot = [max(r["d_sq"], 1e-16) for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.plot(Ns_plot, ds_plot, 'b-', linewidth=1)
    ax.set_xlabel('N'); ax.set_ylabel('d_N')
    ax.set_title('NB Distance vs N')
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.loglog(Ns_plot, ds_plot, 'b-', linewidth=1, label='d_N')
    Ns_arr = np.array(Ns_plot, dtype=float)
    ax.loglog(Ns_arr, A * Ns_arr**(-alpha), 'r--', linewidth=1, label=f'N^(-{alpha:.3f})')
    ax.set_xlabel('N (log)'); ax.set_ylabel('d_N (log)')
    ax.set_title('Power Law Fit')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.semilogy(Ns_plot, dsq_plot, 'g-', linewidth=1)
    ax.set_xlabel('N'); ax.set_ylabel('d_N² (log scale)')
    ax.set_title('Distance Squared (→0 iff RH)')
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Nyman-Beurling: d_N ~ {A:.3f}·N^(-{alpha:.3f})', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "nyman_beurling_200.png"), dpi=150)
    plt.close()
    print(f"\nPlot saved. Results saved.")


if __name__ == "__main__":
    main()

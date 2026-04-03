"""
Gram eigenvalue probe — diagnose rank deficiency onset.

Reuses the 200×200 Gram matrix from v3, loads it if saved, otherwise rebuilds.
Plots: min eigenvalue of G_N vs N, and the eigenvalue spectrum at key N.
Goal: determine if rank deficiency at N~70 is real or precision artifact.
"""

import numpy as np
from scipy.integrate import quad
import json
import os
import time
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR   = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")

# Try loading precomputed G from v3 results json; otherwise rebuild
def frac(x):
    return x - np.floor(x)

def G_entry(j, k, tol=1e-12):
    result, _ = quad(
        lambda t: frac(1.0 / (j * t)) * frac(1.0 / (k * t)),
        1e-10, 1.0, limit=600, epsabs=tol, epsrel=tol
    )
    return result

def b_entry(k, tol=1e-12):
    result, _ = quad(
        lambda t: frac(1.0 / (k * t)),
        1e-10, 1.0, limit=600, epsabs=tol, epsrel=tol
    )
    return result

def build_gram_matrix(N_max):
    G = np.zeros((N_max, N_max))
    b = np.zeros(N_max)
    t0 = time.time()
    total = N_max * (N_max + 1) // 2 + N_max
    done = 0

    for k in range(N_max):
        b[k] = b_entry(k + 1)
        done += 1

    for j in range(N_max):
        for k in range(j, N_max):
            G[j, k] = G_entry(j + 1, k + 1)
            G[k, j] = G[j, k]
            done += 1
            if done % 500 == 0:
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 1
                remaining = (total - done) / rate
                print(f"  [{done}/{total}] {remaining:.0f}s remaining", flush=True)

    print(f"  Done. {time.time()-t0:.1f}s", flush=True)
    return G, b


def main():
    print("GRAM EIGENVALUE PROBE")
    print("=" * 60)

    # Use smaller N for this probe — just need to understand the rank deficiency
    N_MAX = 100

    print(f"Building {N_MAX}×{N_MAX} Gram matrix (higher precision)...")
    G, b = build_gram_matrix(N_MAX)

    print("\nMin/max eigenvalue vs N:")
    print(f"{'N':>5}  {'lambda_min':>14}  {'lambda_max':>14}  {'cond':>12}  {'rank<1e-10':>12}  {'rank<1e-8':>10}  {'rank<1e-6':>10}")
    print("-" * 95)

    records = []
    for N in range(1, N_MAX + 1):
        G_n = G[:N, :N]
        evals = np.linalg.eigvalsh(G_n)
        evals_sorted = np.sort(evals)  # ascending
        lmin = evals_sorted[0]
        lmax = evals_sorted[-1]
        cond = lmax / lmin if lmin > 0 else np.inf
        rank_1e10 = np.sum(evals > 1e-10)
        rank_1e8  = np.sum(evals > 1e-8)
        rank_1e6  = np.sum(evals > 1e-6)

        records.append({
            "N": N, "lambda_min": float(lmin), "lambda_max": float(lmax),
            "cond": float(cond), "rank_1e10": int(rank_1e10),
            "rank_1e8": int(rank_1e8), "rank_1e6": int(rank_1e6)
        })

        if N <= 10 or N % 5 == 0:
            print(f"  {N:>3}  {lmin:>14.4e}  {lmax:>14.4e}  {cond:>12.2e}  {rank_1e10:>12}  {rank_1e8:>10}  {rank_1e6:>10}", flush=True)

    # Find onset of rank deficiency for each tolerance
    for tol_name, key in [("1e-10", "rank_1e10"), ("1e-8", "rank_1e8"), ("1e-6", "rank_1e6")]:
        for r in records:
            if r[key] < r["N"]:
                print(f"\nFirst rank deficiency (tol={tol_name}) at N={r['N']} (rank={r[key]})")
                break

    # Print eigenvalue spectra at N=60, 65, 70, 75
    print("\nEigenvalue spectra at key N:")
    for N in [55, 60, 65, 70, 75, 80]:
        if N > N_MAX:
            continue
        G_n = G[:N, :N]
        evals = np.sort(np.linalg.eigvalsh(G_n))
        print(f"\n  N={N}: smallest 8 eigenvalues:")
        for i, e in enumerate(evals[:8]):
            print(f"    [{i}] {e:.6e}")

    # Save
    with open(os.path.join(RESULTS_DIR, "gram_eigenvalue_probe.json"), "w") as f:
        json.dump({"N_max": N_MAX, "records": records}, f, indent=2)

    # Plot: min eigenvalue vs N (log scale)
    Ns = [r["N"] for r in records]
    lmins = [max(r["lambda_min"], 1e-20) for r in records]  # floor for log plot

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.semilogy(Ns, lmins, 'b-', linewidth=1)
    ax.axhline(1e-10, color='r', linestyle='--', linewidth=0.8, label='tol=1e-10')
    ax.axhline(1e-8,  color='orange', linestyle='--', linewidth=0.8, label='tol=1e-8')
    ax.axhline(1e-6,  color='green', linestyle='--', linewidth=0.8, label='tol=1e-6')
    ax.set_xlabel('N'); ax.set_ylabel('λ_min(G_N)')
    ax.set_title('Minimum Gram Eigenvalue vs N')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1]
    r10 = [r["rank_1e10"] for r in records]
    r8  = [r["rank_1e8"] for r in records]
    r6  = [r["rank_1e6"] for r in records]
    ax.plot(Ns, Ns, 'k--', linewidth=0.8, label='full rank')
    ax.plot(Ns, r10, 'r-', linewidth=1, label='rank (tol=1e-10)')
    ax.plot(Ns, r8,  'orange', linewidth=1, label='rank (tol=1e-8)')
    ax.plot(Ns, r6,  'g-', linewidth=1, label='rank (tol=1e-6)')
    ax.set_xlabel('N'); ax.set_ylabel('Effective rank')
    ax.set_title('Gram Matrix Rank vs N (different tolerances)')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "gram_eigenvalue_probe.png"), dpi=150)
    plt.close()
    print("\nDone. Plot saved.")


if __name__ == "__main__":
    main()

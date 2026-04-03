"""
Nyman-Beurling Distance — v4 (higher precision integration)

Bug found in v3: rank deficiency at N~70 was a precision artifact.
With epsabs=1e-12, limit=600, G_N stays full rank through at least N=100.
The v3 d_N=0 for N>=70 were spurious — integration noise created fake zero eigenvalues.

v4 fix: higher precision quad (epsabs=1e-12, limit=600) + truncated SVD.
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
N_MAX = 200
SVD_TOL = 1e-10


def frac(x):
    return x - np.floor(x)


def G_entry(j, k):
    """∫₀¹ {1/(j·t)} {1/(k·t)} dt — high precision"""
    result, _ = quad(
        lambda t: frac(1.0 / (j * t)) * frac(1.0 / (k * t)),
        1e-10, 1.0, limit=600, epsabs=1e-12, epsrel=1e-12
    )
    return result


def b_entry(k):
    """∫₀¹ {1/(k·t)} dt — high precision"""
    result, _ = quad(
        lambda t: frac(1.0 / (k * t)),
        1e-10, 1.0, limit=600, epsabs=1e-12, epsrel=1e-12
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
                print(f"  [{done}/{total}] (j={j+1},k={k+1}) {remaining:.0f}s remaining", flush=True)

    print(f"  Done. {time.time()-t0:.1f}s", flush=True)
    return G, b


def stable_d_sq(G_n, b_n):
    """
    d²_N = 1 - b^T G^{-1} b via truncated eigendecomposition.
    Returns (d_sq, rank_used, cond_number, lambda_min).
    """
    evals, evecs = np.linalg.eigh(G_n)
    rank = np.sum(evals > SVD_TOL)
    lambda_min = evals[0]

    if rank == 0:
        return 1.0, 0, np.inf, float(evals[0])

    b_proj = evecs.T @ b_n
    kept = evals[evals > SVD_TOL]
    cond = evals[-1] / kept[0] if len(kept) > 0 else np.inf

    d_sq = 1.0
    for i in range(len(evals)):
        if evals[i] > SVD_TOL:
            d_sq -= b_proj[i]**2 / evals[i]

    return max(d_sq, 0.0), rank, cond, float(lambda_min)


def main():
    print("NYMAN-BEURLING DISTANCE v4 — high precision integration, N up to 200")
    print("=" * 70)
    print(f"Integration: epsabs=1e-12, epsrel=1e-12, limit=600")
    print(f"SVD tolerance: {SVD_TOL:.0e}")
    print()

    print(f"Building {N_MAX}×{N_MAX} Gram matrix...")
    G, b = build_gram_matrix(N_MAX)

    print("\nSolving d_N for N=1..200 via truncated SVD...")
    results = []
    prev_d = 1.0
    violations = 0

    for N in range(1, N_MAX + 1):
        G_n = G[:N, :N]
        b_n = b[:N]
        d_sq, rank, cond, lmin = stable_d_sq(G_n, b_n)
        d = np.sqrt(d_sq)

        if d > prev_d + 1e-8:
            violations += 1
        prev_d = d

        results.append({
            "N": N,
            "d_sq": float(d_sq),
            "d": float(d),
            "rank": int(rank),
            "cond": float(cond),
            "lambda_min": float(lmin)
        })

        if N <= 5 or N % 10 == 0:
            rank_str = f"rank={rank}/{N}" if rank < N else f"rank={rank}/{N} (full)"
            print(f"  N={N:>3d}: d_N={d:.8f}  {rank_str}  λ_min={lmin:.3e}  cond={cond:.2e}", flush=True)

    print(f"\nMonotonicity violations: {violations}")

    # Power law fits
    full_rank = [(r["N"], r["d"]) for r in results if r["rank"] == r["N"] and r["d"] > 1e-12]
    if len(full_rank) >= 10:
        Ns_v = np.array([x[0] for x in full_rank], dtype=float)
        ds_v = np.array([x[1] for x in full_rank])
        coeffs = np.polyfit(np.log(Ns_v), np.log(ds_v), 1)
        alpha = -coeffs[0]
        A = np.exp(coeffs[1])
        print(f"\nPower law (all full-rank): d_N ~ {A:.4f} * N^(-{alpha:.4f})")

        third = len(Ns_v) * 2 // 3
        if third > 5:
            coeffs2 = np.polyfit(np.log(Ns_v[third:]), np.log(ds_v[third:]), 1)
            alpha2 = -coeffs2[0]
            A2 = np.exp(coeffs2[1])
            print(f"Power law (asymptotic, N>{Ns_v[third]:.0f}): d_N ~ {A2:.4f} * N^(-{alpha2:.4f})")

    # All points (including rank-deficient)
    all_pts = [(r["N"], r["d"]) for r in results if r["d"] > 1e-12]
    if len(all_pts) >= 10:
        Ns_a = np.array([x[0] for x in all_pts], dtype=float)
        ds_a = np.array([x[1] for x in all_pts])
        coeffs_a = np.polyfit(np.log(Ns_a), np.log(ds_a), 1)
        alpha_a = -coeffs_a[0]
        A_a = np.exp(coeffs_a[1])
        print(f"Power law (all nonzero d_N):  d_N ~ {A_a:.4f} * N^(-{alpha_a:.4f})")

    # Summary
    print("\nSummary:")
    key = [1, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200]
    for r in results:
        if r["N"] in key:
            dep = "" if r["rank"] == r["N"] else f" [rank {r['rank']}/{r['N']}]"
            print(f"  N={r['N']:>3d}: d_N={r['d']:.8f}  λ_min={r['lambda_min']:.3e}{dep}")

    # Save
    with open(os.path.join(RESULTS_DIR, "nyman_beurling_v4.json"), "w") as f:
        json.dump({"N_max": N_MAX, "svd_tol": SVD_TOL, "results": results}, f, indent=2)

    # Plot
    Ns_all = [r["N"] for r in results]
    ds_all = [r["d"] for r in results]
    lmins  = [max(r["lambda_min"], 1e-20) for r in results]
    fr_Ns  = [r["N"] for r in results if r["rank"] == r["N"]]
    fr_ds  = [r["d"] for r in results if r["rank"] == r["N"]]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.plot(Ns_all, ds_all, 'b-', linewidth=0.8, alpha=0.5, label='all N')
    ax.plot(fr_Ns, fr_ds, 'bo', markersize=3, label='full rank')
    ax.set_xlabel('N'); ax.set_ylabel('d_N')
    ax.set_title('NB Distance v4 (high precision)')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1]
    if len(fr_Ns) > 5:
        ax.loglog(fr_Ns, fr_ds, 'bo-', markersize=3, linewidth=0.8)
        ax.set_xlabel('N (log)'); ax.set_ylabel('d_N (log)')
        ax.set_title('Power Law (full rank only)')
        ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.semilogy(Ns_all, lmins, 'g-', linewidth=0.8)
    ax.axhline(SVD_TOL, color='r', linestyle='--', linewidth=0.8, label=f'SVD_TOL={SVD_TOL:.0e}')
    ax.set_xlabel('N'); ax.set_ylabel('λ_min(G_N)')
    ax.set_title('Min Gram Eigenvalue vs N')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "nyman_beurling_v4.png"), dpi=150)
    plt.close()
    print("\nDone.")


if __name__ == "__main__":
    main()

"""
Nyman-Beurling Distance — v5 (breakpoints fix + analytic b_k)

Fixes v4 negative eigenvalue problem:
1. b_k = (1 + log k - γ)/k  — exact analytic formula, no integration
2. G_entry(j,k) uses explicit breakpoints at t=1/(jm) and t=1/(km)
   — these are exactly where the fractional-part functions jump
   — eliminates adaptive subdivision failures for large j,k

Validated: b_1 = 1-γ ≈ 0.4228, G_11 = log(2π)-1-γ ≈ 0.2607.
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
GAMMA = 0.5772156649015328606  # Euler-Mascheroni constant


def frac(x):
    return x - np.floor(x)


def b_entry_analytic(k):
    """Exact: b_k = (1 + log k - γ) / k"""
    return (1.0 + np.log(k) - GAMMA) / k


def G_entry_with_breakpoints(j, k, eps=1e-10):
    """
    ∫₀¹ {1/(j·t)} {1/(k·t)} dt with breakpoints at t=1/(j*m) and t=1/(k*m).

    Breakpoints are where the fractional-part functions jump — providing them
    to quad eliminates the adaptive subdivision failure for large j, k.
    """
    # Collect breakpoints: t = 1/(j*m) for m=1..j, t=1/(k*m) for m=1..k
    # (1/(j*(j+1)) ≈ 1/(j²) — below this, both functions oscillate but
    #  with very small amplitude; use eps as lower cutoff)

    bpts = set()
    # Primary breakpoints (most important — first few periods)
    max_m = max(j, k) + 1
    for m in range(1, max_m + 1):
        p = 1.0 / (j * m)
        if p > eps:
            bpts.add(p)
        p = 1.0 / (k * m)
        if p > eps:
            bpts.add(p)

    bpts = sorted(b for b in bpts if eps < b < 1.0)

    # Scipy quad needs points STRICTLY inside the integration interval
    # and the list must be sorted. Limit must be >= len(points)+1.
    limit = max(600, len(bpts) + 100)

    integrand = lambda t: frac(1.0 / (j * t)) * frac(1.0 / (k * t))

    if len(bpts) == 0:
        result, _ = quad(integrand, eps, 1.0, limit=limit, epsabs=1e-13, epsrel=1e-13)
    else:
        result, _ = quad(integrand, eps, 1.0, limit=limit, epsabs=1e-13, epsrel=1e-13,
                         points=bpts)
    return result


def build_gram_matrix(N_max):
    G = np.zeros((N_max, N_max))
    b = np.zeros(N_max)
    t0 = time.time()
    total = N_max * (N_max + 1) // 2 + N_max
    done = 0

    # b vector: analytic formula, instant
    for k in range(N_max):
        b[k] = b_entry_analytic(k + 1)
        done += 1

    print(f"  b vector: analytic, done instantly.", flush=True)

    for j in range(N_max):
        for k in range(j, N_max):
            G[j, k] = G_entry_with_breakpoints(j + 1, k + 1)
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
    evals, evecs = np.linalg.eigh(G_n)
    rank = np.sum(evals > SVD_TOL)
    lambda_min = float(evals[0])

    if rank == 0:
        return 1.0, 0, np.inf, lambda_min

    b_proj = evecs.T @ b_n
    kept = evals[evals > SVD_TOL]
    cond = float(evals[-1] / kept[0])

    d_sq = 1.0
    for i in range(len(evals)):
        if evals[i] > SVD_TOL:
            d_sq -= b_proj[i]**2 / evals[i]

    return max(d_sq, 0.0), int(rank), cond, lambda_min


def main():
    print("NYMAN-BEURLING DISTANCE v5 — breakpoints fix + analytic b_k")
    print("=" * 70)
    print(f"b_k: analytic formula (1 + log k - γ)/k")
    print(f"G_jk: quad with explicit breakpoints at t=1/(jm), t=1/(km)")
    print(f"SVD tolerance: {SVD_TOL:.0e}")
    print()

    # Validate analytic formulas
    b1_analytic = b_entry_analytic(1)
    G11_expected = np.log(2 * np.pi) - 1 - GAMMA  # 0.2607
    print(f"Validation: b_1 = {b1_analytic:.6f} (expect {1-GAMMA:.6f})")
    print(f"Validation: G_11 expected = {G11_expected:.6f}")
    print()

    print(f"Building {N_MAX}×{N_MAX} Gram matrix...")
    G, b = build_gram_matrix(N_MAX)

    print(f"\nG_11 computed = {G[0,0]:.6f} (expect {G11_expected:.6f})")
    print()

    print("Solving d_N for N=1..200 via truncated SVD...")
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
            "N": N, "d_sq": float(d_sq), "d": float(d),
            "rank": rank, "cond": cond, "lambda_min": lmin
        })

        if N <= 5 or N % 10 == 0:
            full = "(full)" if rank == N else f"[DEFICIENT]"
            print(f"  N={N:>3d}: d_N={d:.8f}  rank={rank}/{N}{full}  λ_min={lmin:.3e}  cond={cond:.2e}", flush=True)

    print(f"\nMonotonicity violations: {violations}")

    # Power law fits
    full_rank_pts = [(r["N"], r["d"]) for r in results if r["rank"] == r["N"] and r["d"] > 1e-12]
    if len(full_rank_pts) >= 10:
        Ns = np.array([x[0] for x in full_rank_pts], dtype=float)
        ds = np.array([x[1] for x in full_rank_pts])
        coeffs = np.polyfit(np.log(Ns), np.log(ds), 1)
        alpha = -coeffs[0]; A = np.exp(coeffs[1])
        print(f"\nPower law (all full-rank):  d_N ~ {A:.4f} * N^(-{alpha:.4f})")

        third = len(Ns) * 2 // 3
        if third > 5:
            c2 = np.polyfit(np.log(Ns[third:]), np.log(ds[third:]), 1)
            a2 = -c2[0]; A2 = np.exp(c2[1])
            print(f"Power law (N>{Ns[third]:.0f}, asymptotic): d_N ~ {A2:.4f} * N^(-{a2:.4f})")

    # Summary
    print("\nSummary:")
    for r in results:
        if r["N"] in [1, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200]:
            dep = "" if r["rank"] == r["N"] else f" [rank {r['rank']}/{r['N']}]"
            print(f"  N={r['N']:>3d}: d_N={r['d']:.8f}  λ_min={r['lambda_min']:.3e}{dep}")

    # Save
    with open(os.path.join(RESULTS_DIR, "nyman_beurling_v5.json"), "w") as f:
        json.dump({"N_max": N_MAX, "svd_tol": SVD_TOL, "results": results}, f, indent=2)

    # Plot
    Ns_all = [r["N"] for r in results]
    ds_all = [r["d"] for r in results]
    lmins  = [max(r["lambda_min"], 1e-20) for r in results]
    fr_N   = [r["N"] for r in results if r["rank"] == r["N"]]
    fr_d   = [r["d"] for r in results if r["rank"] == r["N"]]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.plot(Ns_all, ds_all, 'b-', linewidth=0.8, alpha=0.5)
    ax.plot(fr_N, fr_d, 'bo', markersize=3, label='full rank')
    ax.set_xlabel('N'); ax.set_ylabel('d_N')
    ax.set_title('NB Distance v5 (breakpoints)')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1]
    if len(fr_N) > 5:
        ax.loglog(fr_N, fr_d, 'bo-', markersize=3, linewidth=0.8)
        ax.set_xlabel('N (log)'); ax.set_ylabel('d_N (log)')
        ax.set_title('Power Law (full rank)')
        ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.semilogy(Ns_all, lmins, 'g-', linewidth=0.8)
    ax.axhline(SVD_TOL, color='r', linestyle='--', linewidth=0.8, label='SVD_TOL')
    ax.axhline(0, color='k', linestyle='-', linewidth=0.5)
    ax.set_xlabel('N'); ax.set_ylabel('λ_min(G_N)')
    ax.set_title('Min Gram Eigenvalue vs N')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "nyman_beurling_v5.png"), dpi=150)
    plt.close()
    print("\nDone.")


if __name__ == "__main__":
    main()

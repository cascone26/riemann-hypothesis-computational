"""
Nyman-Beurling Distance — v6 (extended N=300, saves G matrix)

Same method as v5 (analytic b_k, breakpoints for G_{jk}).
Extends to N=300 to get a cleaner asymptotic power-law estimate.
Saves G and b to disk so Báez-Duarte coefficients can be extracted later.

Key question: is α = 0.308 converging toward 1/3, or something else?
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
N_MAX = 300
SVD_TOL = 1e-10
GAMMA = 0.5772156649015328606


def frac(x):
    return x - np.floor(x)


def b_entry_analytic(k):
    return (1.0 + np.log(k) - GAMMA) / k


def G_entry_with_breakpoints(j, k, eps=1e-10):
    max_m = max(j, k) + 1
    bpts = set()
    for m in range(1, max_m + 1):
        for n in (j, k):
            p = 1.0 / (n * m)
            if eps < p < 1.0:
                bpts.add(p)
    bpts = sorted(bpts)
    limit = max(600, len(bpts) + 100)
    integrand = lambda t: frac(1.0 / (j * t)) * frac(1.0 / (k * t))
    if bpts:
        result, _ = quad(integrand, eps, 1.0, limit=limit,
                         epsabs=1e-13, epsrel=1e-13, points=bpts)
    else:
        result, _ = quad(integrand, eps, 1.0, limit=limit,
                         epsabs=1e-13, epsrel=1e-13)
    return result


def build_gram_matrix(N_max):
    G = np.zeros((N_max, N_max))
    b = np.zeros(N_max)
    t0 = time.time()
    total = N_max * (N_max + 1) // 2 + N_max
    done = 0

    for k in range(N_max):
        b[k] = b_entry_analytic(k + 1)
        done += 1

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
    lmin = float(evals[0])
    if rank == 0:
        return 1.0, 0, np.inf, lmin, None
    kept = evals[evals > SVD_TOL]
    cond = float(evals[-1] / kept[0])
    b_proj = evecs.T @ b_n
    d_sq = 1.0
    for i in range(len(evals)):
        if evals[i] > SVD_TOL:
            d_sq -= b_proj[i]**2 / evals[i]
    # Compute coefficients c = G^{-1} b (only for full-rank cases)
    if rank == len(evals):
        # Full rank: compute c = G^{-1} b directly
        c = evecs @ (b_proj / evals)
    else:
        c = None
    return max(d_sq, 0.0), int(rank), cond, lmin, c


def main():
    print("NYMAN-BEURLING DISTANCE v6 — N up to 300")
    print("=" * 70)
    print(f"b_k: analytic (1+log k-γ)/k | G_jk: breakpoints | SVD_TOL={SVD_TOL:.0e}")
    print()

    print(f"Building {N_MAX}×{N_MAX} Gram matrix...")
    G, b = build_gram_matrix(N_MAX)

    # Save G and b for later coefficient analysis
    np.save(os.path.join(RESULTS_DIR, "gram_300.npy"), G)
    np.save(os.path.join(RESULTS_DIR, "bvec_300.npy"), b)
    print(f"G and b saved to results/")

    print(f"G_11 = {G[0,0]:.6f} (analytic: {np.log(2*np.pi)-1-GAMMA:.6f})")
    print()

    print("Solving d_N for N=1..300...")
    results = []
    prev_d = 1.0
    violations = 0
    # Save coefficients at these N values
    coeff_Ns = [50, 100, 150, 200, 250, 300]
    coefficients = {}

    for N in range(1, N_MAX + 1):
        G_n = G[:N, :N]
        b_n = b[:N]
        d_sq, rank, cond, lmin, c = stable_d_sq(G_n, b_n)
        d = np.sqrt(d_sq)

        if d > prev_d + 1e-8:
            violations += 1
        prev_d = d

        results.append({
            "N": N, "d_sq": float(d_sq), "d": float(d),
            "rank": rank, "cond": cond, "lambda_min": lmin
        })

        if N in coeff_Ns and c is not None:
            coefficients[N] = c.tolist()

        if N <= 5 or N % 20 == 0:
            full = "(full)" if rank == N else "[DEFICIENT]"
            print(f"  N={N:>3d}: d_N={d:.8f}  rank={rank}/{N}{full}  λ_min={lmin:.3e}  cond={cond:.2e}", flush=True)

    print(f"\nMonotonicity violations: {violations}")

    # Power law fits at different N ranges
    fr = [(r["N"], r["d"]) for r in results if r["rank"] == r["N"] and r["d"] > 1e-12]
    Ns = np.array([x[0] for x in fr], dtype=float)
    ds = np.array([x[1] for x in fr])

    print("\nPower law fits (d_N ~ A * N^{-α}):")
    for N_lo in [10, 50, 100, 150, 200]:
        mask = Ns >= N_lo
        if mask.sum() < 5:
            continue
        c_fit = np.polyfit(np.log(Ns[mask]), np.log(ds[mask]), 1)
        alpha = -c_fit[0]; A = np.exp(c_fit[1])
        print(f"  N={N_lo}..{int(Ns[-1])}: d_N ~ {A:.4f} * N^(-{alpha:.4f})")

    # Summary
    print("\nSummary:")
    key = [1, 5, 10, 25, 50, 100, 150, 200, 250, 300]
    for r in results:
        if r["N"] in key:
            dep = "" if r["rank"] == r["N"] else f" [rank {r['rank']}/{r['N']}]"
            print(f"  N={r['N']:>3d}: d_N={r['d']:.8f}  λ_min={r['lambda_min']:.3e}{dep}")

    # Báez-Duarte coefficient analysis
    print("\nBáez-Duarte coefficients c_k at selected N:")
    for N_c in sorted(coefficients.keys()):
        c = np.array(coefficients[N_c])
        print(f"  N={N_c}: c_1={c[0]:.6f}  c_2={c[1]:.6f}  c_3={c[2]:.6f}  c_max={c.max():.4f}  c_min={c.min():.4f}")
        # Check if c_k ~ μ(k)/k or similar pattern
        ks = np.arange(1, N_c + 1)
        # Fit c_k vs k: look for power law
        pos = c > 0
        if pos.sum() > 5:
            try:
                cf = np.polyfit(np.log(ks[pos]), np.log(c[pos]), 1)
                print(f"    Power law (positive c_k): c_k ~ {np.exp(cf[1]):.4f} * k^({cf[0]:.3f})")
            except:
                pass

    # Save
    with open(os.path.join(RESULTS_DIR, "nyman_beurling_v6.json"), "w") as f:
        json.dump({
            "N_max": N_MAX, "svd_tol": SVD_TOL,
            "results": results,
            "coefficients": {str(k): v for k, v in coefficients.items()}
        }, f, indent=2)

    # Plot
    Ns_all = [r["N"] for r in results]
    ds_all = [r["d"] for r in results]
    lmins  = [max(r["lambda_min"], 1e-20) for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.plot(Ns_all, ds_all, 'b-', linewidth=1)
    ax.set_xlabel('N'); ax.set_ylabel('d_N')
    ax.set_title('NB Distance v6 (N=1..300)')
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.loglog(Ns_all, ds_all, 'b-', linewidth=1)
    # Power law overlay
    c_fit = np.polyfit(np.log(Ns[Ns >= 100]), np.log(ds[Ns >= 100]), 1)
    alpha_asym = -c_fit[0]; A_asym = np.exp(c_fit[1])
    Ns_fit = np.array([100, 300], dtype=float)
    ax.loglog(Ns_fit, A_asym * Ns_fit**(-alpha_asym), 'r--', linewidth=1,
              label=f'N^(-{alpha_asym:.3f})')
    ax.set_xlabel('N (log)'); ax.set_ylabel('d_N (log)')
    ax.set_title('Power Law Fit')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.semilogy(Ns_all, lmins, 'g-', linewidth=1)
    ax.axhline(SVD_TOL, color='r', linestyle='--', linewidth=0.8, label='SVD_TOL')
    ax.set_xlabel('N'); ax.set_ylabel('λ_min(G_N)')
    ax.set_title('Min Gram Eigenvalue vs N')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "nyman_beurling_v6.png"), dpi=150)
    plt.close()
    print("\nDone.")


if __name__ == "__main__":
    main()

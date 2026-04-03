"""
Nyman-Beurling Distance — v8 (N=2000, same analytic G_{jk} as v7)

Extends v7 to N=2000 to clarify asymptotic rate.
Key question: does alpha (local power law exponent) stabilize, decrease to 0,
or plateau? N=1000 shows alpha ~ 0.02-0.10 in local windows.

Runtime estimate: ~5 min (O(N^3) scaling, 40s at N=1000).
"""

import numpy as np
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
N_MAX  = 2000
SVD_TOL = 1e-10
GAMMA = 0.5772156649015328606


def b_entry_analytic(k):
    return (1.0 + np.log(k) - GAMMA) / k


def G_entry_analytic(j, k):
    M_cutoff = max(k + 1, 100_000 // j + 1)
    m_arr = np.arange(0, M_cutoff + 1, dtype=np.int64)
    n_lo = (j * m_arr) // k
    n_hi = (j * m_arr + j - 1) // k
    max_delta = int((n_hi - n_lo).max())
    total = 0.0
    for delta in range(max_delta + 1):
        n_arr = n_lo + delta
        valid = n_arr <= n_hi
        if not valid.any():
            break
        m_v = m_arr[valid]; n_v = n_arr[valid]
        a = np.maximum(j * m_v, k * n_v).astype(np.float64)
        a = np.maximum(a, 1.0)
        b = np.minimum(j * (m_v + 1), k * (n_v + 1)).astype(np.float64)
        good = b > a
        if not good.any():
            continue
        a, b = a[good], b[good]
        m_f = m_v[good].astype(np.float64); n_f = n_v[good].astype(np.float64)
        contrib = ((b - a) / (j * k)
                   - (n_f / j + m_f / k) * np.log(b / a)
                   + m_f * n_f * (1.0 / a - 1.0 / b))
        total += float(contrib.sum())
    return total


def build_gram_matrix(N_max, report_every=5000):
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
            G[j, k] = G_entry_analytic(j + 1, k + 1)
            G[k, j] = G[j, k]
            done += 1
            if done % report_every == 0:
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
    if rank == len(evals):
        c = evecs @ (b_proj / evals)
    else:
        c = None
    return max(d_sq, 0.0), int(rank), cond, lmin, c


def main():
    print("NYMAN-BEURLING DISTANCE v8 — analytic G_{jk}, N up to 2000")
    print("=" * 75)
    print(f"G_{{jk}}: u-domain cell summation | SVD_TOL={SVD_TOL:.0e}")
    print()

    print(f"Building {N_MAX}×{N_MAX} Gram matrix...")
    G, b = build_gram_matrix(N_MAX)

    np.save(os.path.join(RESULTS_DIR, "gram_2000.npy"), G)
    np.save(os.path.join(RESULTS_DIR, "bvec_2000.npy"), b)
    print(f"G and b saved.")
    print()

    print(f"Solving d_N for N=1..{N_MAX}...")
    results = []
    prev_d = 1.0
    violations = 0
    coeff_Ns = {200, 500, 1000, 1500, 2000}
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

        if N <= 5 or N % 100 == 0:
            full = "(full)" if rank == N else "[DEFICIENT]"
            print(f"  N={N:>4d}: d_N={d:.8f}  rank={rank}/{N}{full}  λ_min={lmin:.3e}  cond={cond:.2e}", flush=True)

    print(f"\nMonotonicity violations: {violations}")

    fr = [(r["N"], r["d"]) for r in results if r["rank"] == r["N"] and r["d"] > 1e-12]
    Ns = np.array([x[0] for x in fr], dtype=float)
    ds = np.array([x[1] for x in fr])

    print("\nPower law fits (d_N ~ A * N^{-α}):")
    for N_lo in [100, 200, 500, 1000, 1500]:
        mask = Ns >= N_lo
        if mask.sum() < 5:
            continue
        c_fit = np.polyfit(np.log(Ns[mask]), np.log(ds[mask]), 1)
        alpha = -c_fit[0]; A = np.exp(c_fit[1])
        print(f"  N={N_lo:>4d}..{int(Ns[-1])}: d_N ~ {A:.4f} * N^(-{alpha:.6f})  [n_pts={mask.sum()}]")

    print("\nRunning α (windows of 100):")
    for N_lo in range(100, N_MAX - 50, 100):
        mask = (Ns >= N_lo) & (Ns < N_lo + 100)
        if mask.sum() < 20:
            continue
        c_fit = np.polyfit(np.log(Ns[mask]), np.log(ds[mask]), 1)
        alpha = -c_fit[0]
        print(f"  N={N_lo:>4d}..{N_lo+99}: α = {alpha:.6f}  d_mid={ds[mask][len(ds[mask])//2]:.6f}")

    print("\nSummary:")
    key_Ns = [1, 5, 10, 25, 50, 100, 200, 300, 500, 750, 1000, 1250, 1500, 1750, 2000]
    for r in results:
        if r["N"] in key_Ns:
            dep = "" if r["rank"] == r["N"] else f" [rank {r['rank']}/{r['N']}]"
            print(f"  N={r['N']:>4d}: d_N={r['d']:.8f}  λ_min={r['lambda_min']:.3e}{dep}")

    with open(os.path.join(RESULTS_DIR, "nyman_beurling_v8.json"), "w") as f:
        json.dump({
            "N_max": N_MAX, "svd_tol": SVD_TOL,
            "results": results,
            "coefficients": {str(k): v for k, v in coefficients.items()}
        }, f, indent=2)

    Ns_all = [r["N"] for r in results]
    ds_all = [r["d"] for r in results]
    lmins  = [max(r["lambda_min"], 1e-20) for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    ax = axes[0]
    ax.plot(Ns_all, ds_all, 'b-', linewidth=0.5)
    ax.set_xlabel('N'); ax.set_ylabel('d_N')
    ax.set_title('NB Distance v8 (N=1..2000, analytic)')
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.loglog(Ns_all, ds_all, 'b-', linewidth=0.5)
    Nref = np.array([200.0, 2000.0])
    mask200 = Ns >= 200
    if mask200.sum() > 5:
        cf = np.polyfit(np.log(Ns[mask200]), np.log(ds[mask200]), 1)
        A_asym = np.exp(cf[1]); alpha_asym = -cf[0]
        ax.loglog(Nref, A_asym * Nref**(-alpha_asym), 'r--', linewidth=1.2,
                  label=f'fit: N^(-{alpha_asym:.4f})')
    ax.loglog(Nref, 0.5 * Nref**(-1/3), 'g:', linewidth=1, label='N^(-1/3) ref')
    ax.set_xlabel('N (log)'); ax.set_ylabel('d_N (log)')
    ax.set_title('Power Law Fit')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.semilogy(Ns_all, lmins, 'g-', linewidth=0.5)
    ax.axhline(SVD_TOL, color='r', linestyle='--', linewidth=0.8, label='SVD_TOL')
    ax.set_xlabel('N'); ax.set_ylabel('λ_min(G_N)')
    ax.set_title('Min Gram Eigenvalue vs N')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "nyman_beurling_v8.png"), dpi=150)
    plt.close()
    print("\nDone.")


if __name__ == "__main__":
    main()

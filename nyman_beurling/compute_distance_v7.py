"""
Nyman-Beurling Distance — v7 (fully analytic G_{jk}, u-domain summation)

Key improvement over v6: replaces scipy.integrate.quad with exact u-domain
cell summation (u = 1/t substitution).

G_{jk} = ∫₁^∞ {u/j}{u/k}/u² du

Summed over cells where floor(u/j)=m and floor(u/k)=n are both constant.
For j ≤ k, each m-cell has at most 2 n-values (delta=0 or 1).

Cell [a,b]: a=max(j*m, k*n), b=min(j*(m+1), k*(n+1))
Cell integral (exact):
  (b-a)/(j*k) - (n/j + m/k)*log(b/a) + m*n*(1/a - 1/b)

n-range via integer arithmetic:
  n_lo = (j*m)//k,  n_hi = (j*m + j - 1)//k

Tail from m > M_cutoff: bounded by 1/(4*j*M_cutoff) ~ 3e-6 absolute.
Speed: ~10-50x faster than v6, enables N=1000 in a few minutes.
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
N_MAX  = 1000
SVD_TOL = 1e-10
GAMMA = 0.5772156649015328606


def frac(x):
    return x - np.floor(x)


def b_entry_analytic(k):
    return (1.0 + np.log(k) - GAMMA) / k


def G_entry_analytic(j, k):
    """
    G_{jk} = ∫₁^∞ {u/j}{u/k}/u² du  (via u=1/t substitution)

    Sums over cells (m,n) where floor(u/j)=m, floor(u/k)=n are constant.
    Assumes j <= k.  For j<=k, each m-cell has at most 2 n-values (delta=0 or 1).

    m starts from 0: when j>1, u ∈ [1,j) has m=0 and contributes ∫₁^j {u/j}{u/k}/u² du.
    Lower cell boundary is clamped to 1 since the integral starts at u=1.

    n-range via integer arithmetic:  n_lo = (j*m)//k,  n_hi = (j*m+j-1)//k
    Cell integral (exact): (b-a)/(j*k) - (n/j + m/k)*log(b/a) + m*n*(1/a - 1/b)
    Tail from m > M_cutoff bounded by 1/(4*j*M_cutoff) in absolute value.
    """
    M_cutoff = max(k + 1, 100_000 // j + 1)

    m_arr = np.arange(0, M_cutoff + 1, dtype=np.int64)   # include m=0
    n_lo = (j * m_arr) // k
    n_hi = (j * m_arr + j - 1) // k
    max_delta = int((n_hi - n_lo).max())   # 0 or 1 when j <= k

    total = 0.0
    for delta in range(max_delta + 1):
        n_arr = n_lo + delta
        valid = n_arr <= n_hi
        if not valid.any():
            break

        m_v = m_arr[valid]
        n_v = n_arr[valid]

        # Clamp lower bound to 1: integral starts at u=1, not u=j*m
        a = np.maximum(j * m_v, k * n_v).astype(np.float64)
        a = np.maximum(a, 1.0)
        b = np.minimum(j * (m_v + 1), k * (n_v + 1)).astype(np.float64)

        good = b > a
        if not good.any():
            continue
        a, b = a[good], b[good]
        m_f = m_v[good].astype(np.float64)
        n_f = n_v[good].astype(np.float64)

        contrib = ((b - a) / (j * k)
                   - (n_f / j + m_f / k) * np.log(b / a)
                   + m_f * n_f * (1.0 / a - 1.0 / b))
        total += float(contrib.sum())

    return total


def build_gram_matrix(N_max, report_every=1000):
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
    print("NYMAN-BEURLING DISTANCE v7 — analytic G_{jk}, no scipy, N up to 1000")
    print("=" * 75)
    print(f"G_{'{jk}'}: fully analytic (closed-form sub-intervals, same breakpoints as v5/v6)")
    print(f"b_k: analytic (1+log k-γ)/k | SVD_TOL={SVD_TOL:.0e}")
    print()

    # Validate against known values
    G11 = G_entry_analytic(1, 1)
    G11_exact = np.log(2 * np.pi) - 1 - GAMMA
    print(f"Validation: G_11 = {G11:.8f} (analytic = {G11_exact:.8f}, err = {abs(G11-G11_exact):.2e})")
    print(f"Validation: b_1  = {b_entry_analytic(1):.8f} (expect = {1-GAMMA:.8f})")
    print()

    print(f"Building {N_MAX}×{N_MAX} Gram matrix...")
    G, b = build_gram_matrix(N_MAX)

    # Save matrices
    np.save(os.path.join(RESULTS_DIR, "gram_1000.npy"), G)
    np.save(os.path.join(RESULTS_DIR, "bvec_1000.npy"), b)
    print(f"G and b saved to results/")
    print()

    print("Solving d_N for N=1..1000...")
    results = []
    prev_d = 1.0
    violations = 0
    coeff_Ns = [100, 200, 300, 500, 750, 1000]
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

        if N <= 5 or N % 50 == 0:
            full = "(full)" if rank == N else "[DEFICIENT]"
            print(f"  N={N:>4d}: d_N={d:.8f}  rank={rank}/{N}{full}  λ_min={lmin:.3e}  cond={cond:.2e}", flush=True)

    print(f"\nMonotonicity violations: {violations}")

    # Power law fits at multiple ranges
    fr = [(r["N"], r["d"]) for r in results if r["rank"] == r["N"] and r["d"] > 1e-12]
    Ns = np.array([x[0] for x in fr], dtype=float)
    ds = np.array([x[1] for x in fr])

    print("\nPower law fits (d_N ~ A * N^{-α}):")
    for N_lo in [10, 50, 100, 200, 300, 500, 700]:
        mask = Ns >= N_lo
        if mask.sum() < 5:
            continue
        c_fit = np.polyfit(np.log(Ns[mask]), np.log(ds[mask]), 1)
        alpha = -c_fit[0]; A = np.exp(c_fit[1])
        print(f"  N={N_lo:>4d}..{int(Ns[-1])}: d_N ~ {A:.4f} * N^(-{alpha:.6f})  [n_pts={mask.sum()}]")

    # Running α estimate
    print("\nRunning α estimate (last 100 points):")
    for N_lo in range(100, N_MAX - 50, 50):
        mask = (Ns >= N_lo) & (Ns < N_lo + 100)
        if mask.sum() < 10:
            continue
        c_fit = np.polyfit(np.log(Ns[mask]), np.log(ds[mask]), 1)
        alpha = -c_fit[0]
        print(f"  N={N_lo:>4d}..{N_lo+99}: α = {alpha:.6f}")

    # Summary table
    print("\nSummary:")
    key_Ns = [1, 5, 10, 25, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    for r in results:
        if r["N"] in key_Ns:
            dep = "" if r["rank"] == r["N"] else f" [rank {r['rank']}/{r['N']}]"
            print(f"  N={r['N']:>4d}: d_N={r['d']:.8f}  λ_min={r['lambda_min']:.3e}{dep}")

    # Báez-Duarte coefficient analysis
    if coefficients:
        print("\nBáez-Duarte coefficients c_k at selected N:")
        for N_c in sorted(coefficients.keys()):
            c = np.array(coefficients[N_c])
            print(f"  N={N_c}: c_1={c[0]:.6f}  c_2={c[1]:.6f}  c_3={c[2]:.6f}  |c_max|={np.abs(c).max():.4f}")
            ks = np.arange(1, N_c + 1)
            pos = c > 0
            if pos.sum() > 5:
                try:
                    cf = np.polyfit(np.log(ks[pos]), np.log(c[pos]), 1)
                    print(f"    Power law (c_k > 0): c_k ~ {np.exp(cf[1]):.4f} * k^({cf[0]:.4f})")
                except:
                    pass

    # Save results
    with open(os.path.join(RESULTS_DIR, "nyman_beurling_v7.json"), "w") as f:
        json.dump({
            "N_max": N_MAX, "svd_tol": SVD_TOL,
            "results": results,
            "coefficients": {str(k): v for k, v in coefficients.items()}
        }, f, indent=2)

    # Plots
    Ns_all = [r["N"] for r in results]
    ds_all = [r["d"] for r in results]
    lmins  = [max(r["lambda_min"], 1e-20) for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.plot(Ns_all, ds_all, 'b-', linewidth=0.8)
    ax.set_xlabel('N'); ax.set_ylabel('d_N')
    ax.set_title('NB Distance v7 (N=1..1000, analytic)')
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.loglog(Ns_all, ds_all, 'b-', linewidth=0.8)
    # Overlay α=1/3 reference
    Nref = np.array([100.0, 1000.0])
    # Fit from N>=200
    mask200 = Ns >= 200
    if mask200.sum() > 5:
        cf = np.polyfit(np.log(Ns[mask200]), np.log(ds[mask200]), 1)
        A_asym = np.exp(cf[1]); alpha_asym = -cf[0]
        ax.loglog(Nref, A_asym * Nref**(-alpha_asym), 'r--', linewidth=1.2,
                  label=f'fit: N^(-{alpha_asym:.4f})')
    ax.loglog(Nref, 0.5 * Nref**(-1/3), 'g:', linewidth=1, label='N^(-1/3) ref')
    ax.set_xlabel('N (log)'); ax.set_ylabel('d_N (log)')
    ax.set_title('Power Law Fit — is α = 1/3?')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.semilogy(Ns_all, lmins, 'g-', linewidth=0.8)
    ax.axhline(SVD_TOL, color='r', linestyle='--', linewidth=0.8, label='SVD_TOL')
    ax.set_xlabel('N'); ax.set_ylabel('λ_min(G_N)')
    ax.set_title('Min Gram Eigenvalue vs N')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "nyman_beurling_v7.png"), dpi=150)
    plt.close()
    print("\nDone.")


if __name__ == "__main__":
    main()

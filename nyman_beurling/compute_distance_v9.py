"""
Nyman-Beurling Distance — v9 (sparse N=5000, key N values only)

Same analytic G_{jk} as v7/v8.
Only solves d_N at a sparse set of N values to avoid the O(N^4) solve bottleneck.
Targets N up to 5000. Matrix build: ~40 min. Sparse solve: minutes.
"""

import numpy as np
import json
import os
import time
import warnings

warnings.filterwarnings('ignore')

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
N_MAX  = 5000
SVD_TOL = 1e-10
GAMMA = 0.5772156649015328606

SPARSE_NS = (list(range(1, 20)) +
             list(range(20, 100, 5)) +
             list(range(100, 500, 10)) +
             list(range(500, 1000, 25)) +
             list(range(1000, 2000, 50)) +
             list(range(2000, 5001, 100)))
SPARSE_NS = sorted(set(SPARSE_NS))


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


def build_gram_matrix(N_max, report_every=10000):
    G = np.zeros((N_max, N_max))
    b = np.zeros(N_max)
    t0 = time.time()
    total = N_max * (N_max + 1) // 2
    done = 0

    for k in range(N_max):
        b[k] = b_entry_analytic(k + 1)

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
        return 1.0, 0, np.inf, lmin
    kept = evals[evals > SVD_TOL]
    cond = float(evals[-1] / kept[0])
    b_proj = evecs.T @ b_n
    d_sq = 1.0
    for i in range(len(evals)):
        if evals[i] > SVD_TOL:
            d_sq -= b_proj[i]**2 / evals[i]
    return max(d_sq, 0.0), int(rank), cond, lmin


def main():
    print(f"NYMAN-BEURLING DISTANCE v9 — sparse, N up to {N_MAX}")
    print("=" * 75)
    print(f"Solving at {len(SPARSE_NS)} sparse N values: {SPARSE_NS[:5]}...{SPARSE_NS[-5:]}")
    print()

    print(f"Building {N_MAX}×{N_MAX} Gram matrix...")
    G, b = build_gram_matrix(N_MAX)

    np.save(os.path.join(RESULTS_DIR, "gram_5000.npy"), G)
    np.save(os.path.join(RESULTS_DIR, "bvec_5000.npy"), b)
    print(f"Saved gram_5000.npy and bvec_5000.npy")
    print()

    print("Solving d_N at sparse N values...")
    results = []
    t0 = time.time()
    for i, N in enumerate(SPARSE_NS):
        G_n = G[:N, :N]
        b_n = b[:N]
        d_sq, rank, cond, lmin = stable_d_sq(G_n, b_n)
        d = np.sqrt(d_sq)
        results.append({"N": N, "d_sq": float(d_sq), "d": float(d),
                        "rank": rank, "cond": cond, "lambda_min": lmin})
        full = "(full)" if rank == N else "[DEFICIENT]"
        print(f"  N={N:>5d}: d_N={d:.8f}  {full}  λ_min={lmin:.3e}  cond={cond:.2e}", flush=True)
    print(f"  Solve done in {time.time()-t0:.1f}s")

    # Rate analysis
    Ns_arr = np.array([r["N"] for r in results if r["rank"] == r["N"]], dtype=float)
    ds_arr = np.array([r["d"] for r in results if r["rank"] == r["N"]])

    print("\nPower law fits:")
    for N_lo in [100, 500, 1000, 2000, 3000, 4000]:
        mask = Ns_arr >= N_lo
        if mask.sum() < 3:
            continue
        c = np.polyfit(np.log(Ns_arr[mask]), np.log(ds_arr[mask]), 1)
        print(f"  N={N_lo:>5d}..{int(Ns_arr[-1])}: α={-c[0]:.6f}  A={np.exp(c[1]):.5f}")

    print("\nLog law fits (d ~ A/log(N)^β):")
    from scipy.optimize import curve_fit
    for N_lo in [500, 1000, 2000, 3000]:
        mask = Ns_arr >= N_lo
        if mask.sum() < 3:
            continue
        try:
            popt, _ = curve_fit(lambda N, A, b: A * np.log(N)**(-b),
                                Ns_arr[mask], ds_arr[mask], p0=[0.3, 0.4])
            print(f"  N={N_lo:>5d}..{int(Ns_arr[-1])}: β={popt[1]:.6f}  A={popt[0]:.5f}")
        except:
            pass

    with open(os.path.join(RESULTS_DIR, "nyman_beurling_v9.json"), "w") as f:
        json.dump({"N_max": N_MAX, "svd_tol": SVD_TOL, "results": results}, f, indent=2)
    print("\nDone.")


if __name__ == "__main__":
    main()

"""
NB Worker — standalone Nyman-Beurling d_N computation for a single N.
Usage: python nb_worker.py <N>
Outputs JSON result to stdout.
"""
import sys, time, json
import numpy as np

GAMMA = 0.5772156649015328606
SVD_TOL = 1e-10

def b_entry(k):
    return (1.0 + np.log(k) - GAMMA) / k

def G_entry(j, k):
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
        m_f = m_v[good].astype(np.float64)
        n_f = n_v[good].astype(np.float64)
        contrib = ((b - a) / (j * k)
                   - (n_f / j + m_f / k) * np.log(b / a)
                   + m_f * n_f * (1.0 / a - 1.0 / b))
        total += float(contrib.sum())
    return total

def compute_d(N):
    t0 = time.time()
    G = np.zeros((N, N))
    b = np.array([b_entry(k+1) for k in range(N)])
    total_pairs = N * (N + 1) // 2
    done = 0
    for j in range(N):
        for k in range(j, N):
            v = G_entry(j+1, k+1)
            G[j, k] = v
            G[k, j] = v
            done += 1
            if done % max(1, total_pairs // 20) == 0:
                elapsed = time.time() - t0
                eta = elapsed / done * (total_pairs - done)
                print(f"  [{done}/{total_pairs}] eta={eta:.0f}s", flush=True)

    t1 = time.time()
    evals, evecs = np.linalg.eigh(G)
    b_proj = evecs.T @ b
    d_sq = 1.0
    for i in range(len(evals)):
        if evals[i] > SVD_TOL:
            d_sq -= b_proj[i]**2 / evals[i]
    d_sq = max(0.0, d_sq)
    d = float(np.sqrt(d_sq))
    t2 = time.time()

    lmin = float(evals[0])
    lmax = float(evals[-1])
    cond = lmax / max(lmin, SVD_TOL)
    rank = int(np.sum(evals > SVD_TOL))

    result = {
        "N": N,
        "d_N": d,
        "d_sq": float(d_sq),
        "lambda_min": lmin,
        "lambda_max": lmax,
        "condition": cond,
        "rank": rank,
        "build_time": round(t1 - t0, 2),
        "solve_time": round(t2 - t1, 2),
        "total_time": round(t2 - t0, 2),
    }
    print(json.dumps(result), flush=True)
    return result

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    print(f"Computing d_N for N={N}...", flush=True)
    compute_d(N)

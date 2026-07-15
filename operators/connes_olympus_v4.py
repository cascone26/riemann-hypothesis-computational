#!/usr/bin/env python3
"""
connes_olympus_v4.py — Large-N Connes sweep on Olympus.

Same math as connes_decisive_v3.py (non-circular, verified correct):
   tau[n,m] = (1/2pi) int Omega(r) ehat_n(r) conj(ehat_m(r)) dr   [arch]
            + W_{0,2}[n,m]                                          [boundary]
            - sum_{k<=lam^2} Lambda(k) k^{-1/2} q(n,m,log k)       [prime]
Omega(r) = Re psi(1/4 + i r/2) - log pi   (digamma, ZERO-FREE)

Change vs v3: arch is built in CHUNKS of G_CHUNK grid points so that
memory stays bounded regardless of N.  The sum over chunks is exact.

Sweep: N = 16, 24, 32, 48, 64, 96, 128, 200, 256, 300, 400, 500
Outputs per N: lambda_min, 5 smallest eigenvalues, zero prediction errors.
Streams to JSON after every N so results survive early termination.
"""

import json, time, sys, os
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30

LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
G_CHUNK = 50_000   # grid points per chunk — keeps ~600MB peak per chunk at N=500

OUTFILE = os.path.join(os.path.expanduser("~"), "projects", "results", "connes_olympus_v4.json")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

N_VALUES = [16, 24, 32, 48, 64, 96, 128, 200, 256, 300, 400, 500]

# Reference zeros from mpmath (independent, never from a file)
print("Computing mpmath reference zeros...", flush=True)
REFZEROS = [float(mp.im(mp.zetazero(k))) for k in range(1, 26)]
print(f"  gamma_1..25: {[round(z,4) for z in REFZEROS[:5]]}...", flush=True)


def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab

VM = vm_table(LAMBDA_SQ)
print(f"Von Mangoldt table: {len(VM)} entries up to k={LAMBDA_SQ}", flush=True)


def omega_vec(r):
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI


def ehat_vec(r, n):
    s = TWO_PI * n / L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    sm = np.abs(d) < 1e-12
    out[sm] = np.sqrt(L) + 0j
    dd = d[~sm]
    out[~sm] = (1.0 / np.sqrt(L)) * np.exp(0.5j * dd * L) * 2.0 * np.sin(0.5 * dd * L) / dd
    return out


def q_kernel_vec(n, m, y_arr):
    if n == m:
        return 2.0 * (1.0 - np.abs(y_arr) / L) * np.cos(TWO_PI * n * y_arr / L)
    return (np.sin(TWO_PI * m * y_arr / L) - np.sin(TWO_PI * n * y_arr / L)) / (np.pi * (n - m))


def prime_mat(idx):
    D = len(idx)
    P = np.zeros((D, D))
    log_k = np.array([np.log(k) for k, _ in VM])   # y argument is log(k), not Lambda(k)
    weights = np.array([lk * k**(-0.5) for k, lk in VM])
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            q = np.where(n == m,
                         2.0 * (1.0 - np.abs(log_k) / L) * np.cos(TWO_PI * n * log_k / L),
                         (np.sin(TWO_PI * m * log_k / L) - np.sin(TWO_PI * n * log_k / L)) / (np.pi * (n - m) + 1e-300))
            P[ii, jj] = P[jj, ii] = float(weights @ q)
    return P


def w02_mat(idx):
    D = len(idx)
    sh = np.sinh(L / 4) ** 2
    ns = np.array(idx, dtype=float)
    num = (L**2 - 16 * np.pi**2 * np.outer(ns, ns))
    den = (L**2 + 16 * np.pi**2 * ns**2)
    W = 32 * L * sh * num / (np.outer(den, den))
    return W


def build_arch_chunked(idx, R, G_total):
    D = len(idx)
    arch = np.zeros((D, D))
    n_chunks = max(1, (G_total + G_CHUNK - 1) // G_CHUNK)
    total_w = 0.0
    for c in range(n_chunks):
        g0 = c * G_CHUNK
        g1 = min(g0 + G_CHUNK, G_total)
        grid = np.linspace(-R + 2*R*g0/G_total, -R + 2*R*g1/G_total, g1 - g0, endpoint=(c == n_chunks-1))
        h = 2*R / G_total
        w = np.full(len(grid), h)
        if c == 0:          w[0]  = 0.5 * h
        if c == n_chunks-1: w[-1] = 0.5 * h
        Om = omega_vec(grid)
        E = np.empty((D, len(grid)), dtype=complex)
        for k, n in enumerate(idx):
            E[k] = ehat_vec(grid, n)
        A = E * (w * Om)
        arch += (A @ E.conj().T).real / (2.0 * np.pi)
        total_w += w.sum()
    return arch


def build_tau(N):
    idx = list(range(-N, N + 1))
    D = len(idx)
    R = TWO_PI * N / L + 80.0
    G_total = max(200_000, 3000 * N)
    t0 = time.time()
    print(f"  arch (chunked G={G_total}, R={R:.1f})...", flush=True)
    arch = build_arch_chunked(idx, R, G_total)
    t1 = time.time()
    print(f"    done {t1-t0:.1f}s", flush=True)
    print(f"  prime matrix...", flush=True)
    P = prime_mat(idx)
    t2 = time.time()
    print(f"    done {t2-t1:.1f}s", flush=True)
    W = w02_mat(idx)
    tau = arch + W - P
    tau = (tau + tau.T) / 2
    return tau, idx


def xi_hat_factory(xi, idx):
    poles = np.array([TWO_PI * j / L for j in idx])
    xi = np.array(xi, dtype=float)
    def xh(z):
        d = z - poles
        d[np.abs(d) < 1e-12] = 1e-12
        return (2.0 / np.sqrt(L)) * np.sin(z * L / 2.0) * np.sum(xi / d)
    return xh


def find_zeros(xh, lo, hi, ng):
    xs = np.linspace(lo, hi, ng)
    vs = np.array([xh(x) for x in xs])
    out = []
    for i in range(1, len(xs)):
        if vs[i-1] == 0 or (vs[i-1] < 0) != (vs[i] < 0):
            a, b = xs[i-1], xs[i]
            for _ in range(80):
                mid = 0.5*(a+b)
                if (xh(a)<0) != (xh(mid)<0): b = mid
                else: a = mid
            r = 0.5*(a+b)
            if not out or abs(r - out[-1]) > 1e-4:
                out.append(r)
    return out


def score_zeros(found, nref):
    errs = []
    for gz in REFZEROS[:nref]:
        if found:
            best = min(found, key=lambda z: abs(z - gz))
            errs.append(abs(best - gz))
        else:
            errs.append(float('nan'))
    return errs


def run_N(N):
    print(f"\n{'='*60}", flush=True)
    print(f"N = {N}  (matrix {2*N+1}x{2*N+1})", flush=True)
    t_start = time.time()

    tau, idx = build_tau(N)

    print(f"  eigendecomposition...", flush=True)
    t_eig = time.time()
    ev, evec = np.linalg.eigh(tau)
    print(f"    done {time.time()-t_eig:.1f}s  lam_min={ev[0]:.6e}", flush=True)

    # Smallest eigenvector → secular function zeros
    k0 = int(np.argmin(np.abs(ev)))
    D = 2*N+1
    v = evec[:, k0]
    even = np.array([(v[i] + v[D-1-i]) / 2 for i in range(D)])
    s = even.sum() / np.sqrt(L)
    if abs(s) > 1e-30:
        even /= s
    xh = xi_hat_factory(even, idx)

    # Find zeros of secular function (up to the bandwidth ceiling s_N = 2*pi*N/L)
    bandwidth_ceil = TWO_PI * N / L
    hi_search = min(bandwidth_ceil * 0.9, 300.0)
    ng = max(5000, int(hi_search * 100))
    zeros = find_zeros(xh, 10.0, hi_search, ng)

    # How many reference zeros fall within bandwidth?
    n_ref_in_bw = sum(1 for gz in REFZEROS if gz < bandwidth_ceil)
    errs = score_zeros(zeros, min(n_ref_in_bw, 25))
    tight = sum(1 for e in errs if e is not None and not np.isnan(e) and e < 0.01)
    loose = sum(1 for e in errs if e is not None and not np.isnan(e) and e < 0.5)

    print(f"  Zeros in [{10},{hi_search:.0f}]: found {len(zeros)}, "
          f"tight(<0.01)={tight}, loose(<0.5)={loose}/{n_ref_in_bw}", flush=True)
    print(f"  5 smallest eigvals: {[f'{e:.4e}' for e in ev[:5]]}", flush=True)
    print(f"  Total time: {time.time()-t_start:.1f}s", flush=True)

    return {
        "N": N,
        "dim": 2*N+1,
        "lambda_min": float(ev[0]),
        "5_smallest": ev[:5].tolist(),
        "5_largest": ev[-5:].tolist(),
        "n_zeros_found": len(zeros),
        "n_ref_in_bandwidth": n_ref_in_bw,
        "zeros_tight": tight,
        "zeros_loose": loose,
        "zero_errors_first25": [float(e) if not np.isnan(e) else None for e in errs],
        "bandwidth_ceil": float(bandwidth_ceil),
        "time_s": round(time.time() - t_start, 1)
    }


def main():
    print("CONNES OLYMPUS v4 — Large-N sweep", flush=True)
    print(f"lambda^2={LAMBDA_SQ}, L={L:.4f}", flush=True)
    print(f"N values: {N_VALUES}", flush=True)
    print(f"G_CHUNK={G_CHUNK}", flush=True)
    print(f"Output: {OUTFILE}", flush=True)

    results = []
    out = {
        "description": "Large-N Connes sweep, chunked arch, non-circular v3 construction",
        "lambda_sq": LAMBDA_SQ, "L": L, "G_CHUNK": G_CHUNK,
        "ref_zeros_25": REFZEROS[:25],
        "runs": results
    }

    for N in N_VALUES:
        try:
            r = run_N(N)
            results.append(r)
            with open(OUTFILE, "w") as f:
                json.dump(out, f, indent=2)
            print(f"  Saved -> {OUTFILE}", flush=True)
        except Exception as e:
            print(f"  ERROR at N={N}: {e}", flush=True)
            results.append({"N": N, "error": str(e)})
            with open(OUTFILE, "w") as f:
                json.dump(out, f, indent=2)
            if "memory" in str(e).lower() or "MemoryError" in type(e).__name__:
                print("  OOM — stopping sweep.", flush=True)
                break

    print("\n" + "="*60, flush=True)
    print("FINAL SUMMARY", flush=True)
    for r in results:
        if "error" not in r:
            print(f"  N={r['N']:4d}  lam_min={r['lambda_min']:+.4e}  "
                  f"zeros tight/loose={r['zeros_tight']}/{r['zeros_loose']}  "
                  f"t={r['time_s']}s", flush=True)
    print(f"\nAll results: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()

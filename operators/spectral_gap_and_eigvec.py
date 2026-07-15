#!/usr/bin/env python3
"""
spectral_gap_and_eigvec.py — Spectral gap analysis and near-zero eigenvector structure.

Questions:
  1. Is there a spectral gap between the ~23 near-zero eigenvalues and the rest?
  2. What are the near-zero eigenvectors? Can we identify them with zeta zeros?
  3. How do near-zero eigenvectors project onto prime-power modes?
  4. Is the spectral gap growing, shrinking, or constant with N?
"""

import json, os, time, numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
LAMBDA_SQ = 14
L = float(mp.log(LAMBDA_SQ))
sh = float(np.sinh(L/4)**2)

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/spectral_gap_analysis.json")

REF_ZEROS = [14.134725, 21.022040, 25.010856, 30.424876, 32.935062,
             37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
             52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
             67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
             79.337376, 82.910381, 84.735493, 87.425275, 88.809112]


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


def build_tau(N):
    VM = vm_table(LAMBDA_SQ)
    idx = list(range(-N, N+1))
    D = len(idx)
    R = TWO_PI * N / L + 80.0
    G = max(60_000, 3000 * N)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5*h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI

    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        s = TWO_PI * n / L
        d = s - grid
        sm = np.abs(d) < 1e-12
        E[k, sm] = np.sqrt(L)
        dd = d[~sm]
        E[k, ~sm] = (1/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2*np.sin(0.5*dd*L) / dd

    arch = ((E * (w * Om)) @ E.conj().T).real / (2*np.pi)
    ns = np.array(idx, dtype=float)
    B = 32*L*sh*(L**2 - 16*np.pi**2*np.outer(ns,ns)) / np.outer(L**2+16*np.pi**2*ns**2, L**2+16*np.pi**2*ns**2)
    log_k = np.array([np.log(k) for k, _ in VM])
    wts = np.array([lk * k**(-0.5) for k, lk in VM])
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            q = np.where(n == m,
                2*(1-np.abs(log_k)/L)*np.cos(TWO_PI*n*log_k/L),
                (np.sin(TWO_PI*m*log_k/L)-np.sin(TWO_PI*n*log_k/L))/(np.pi*(n-m)+1e-300))
            P[ii,jj] = P[jj,ii] = float(wts @ q)

    tau = arch + B - P
    return (tau + tau.T)/2, VM


def analyze(N, tau, VM):
    idx = list(range(-N, N+1))
    D = len(idx)

    ev, evec = np.linalg.eigh(tau)
    s_N = TWO_PI * N / L

    # Near-zero cluster
    n_nz = np.sum(np.abs(ev) < 1e-4)
    # Find spectral gap
    positive_ev = ev[ev > 0]
    gap_val = float(positive_ev[n_nz]) if len(positive_ev) > n_nz else 0.0
    # Actually: first eigenvalue above threshold
    above = ev[ev > 1e-4]
    first_positive_above = float(above[0]) if len(above) > 0 else None

    print(f"\n=== N={N}, s_N={s_N:.2f} ===")
    print(f"  Smallest 30 eigenvalues:")
    for i in range(min(30, len(ev))):
        marker = " <-- NEAR-ZERO" if abs(ev[i]) < 1e-4 else ""
        print(f"    [{i:2d}] {ev[i]:12.6e}{marker}")

    print(f"\n  near-zero count (1e-4): {n_nz}")
    print(f"  First eigenvalue above 1e-4: {first_positive_above:.6e}" if first_positive_above else "  (none above 1e-4)")
    print(f"  Spectral gap ratio: {first_positive_above/max(abs(ev[0]),1e-15):.2e}" if first_positive_above else "")

    # Near-zero eigenvectors: project onto n-modes
    nz_vecs = evec[:, :n_nz]  # D×n_nz
    power_per_mode = np.abs(nz_vecs)**2  # D×n_nz
    # For each eigenvector, which n modes carry the most weight?
    print(f"\n  Near-zero eigenvector analysis (top 3 modes per eigvec):")
    for j in range(min(n_nz, 8)):
        pwr = power_per_mode[:, j]
        top3 = np.argsort(pwr)[-3:][::-1]
        mode_str = ", ".join([f"n={idx[i]}({pwr[i]:.3f})" for i in top3])
        print(f"    evec[{j}] lambda={ev[j]:.2e}: {mode_str}")

    # Which zeta zeros match the dominant n-modes?
    # s_n = 2*pi*n/L, match to gamma_k
    print(f"\n  Zeta zero mode matching (s_n = 2π n/L):")
    for k, gamma in enumerate(REF_ZEROS):
        n_match = round(gamma * L / (TWO_PI))
        s_match = TWO_PI * n_match / L
        if abs(n_match) <= N:
            err = abs(s_match - gamma)
            print(f"    gamma_{k+1}={gamma:.3f}: closest n={n_match}, s_n={s_match:.3f}, err={err:.3f}")

    # Compute projection of each near-zero eigvec onto zero-aligned basis
    zero_n_list = [round(g * L / TWO_PI) for g in REF_ZEROS if abs(round(g * L / TWO_PI)) <= N]
    print(f"\n  Projection of near-zero subspace onto zero-aligned modes:")
    zero_proj = 0.0
    for nz_idx in zero_n_list:
        i_idx = nz_idx + N
        if 0 <= i_idx < D:
            col = nz_vecs[i_idx, :]  # component of each near-zero vec in mode n=nz_idx
            zero_proj += np.sum(col**2)
    total_proj = np.sum(nz_vecs**2)
    print(f"    {zero_proj/total_proj*100:.1f}% of near-zero subspace on zero-aligned modes")

    return {
        'N': N,
        's_N': float(s_N),
        'lambda_min': float(ev[0]),
        'lambda_max': float(ev[-1]),
        'n_near_zero_1e4': int(n_nz),
        'first_above_1e4': first_positive_above,
        'spectral_gap_ratio': float(first_positive_above / max(abs(ev[0]), 1e-15)) if first_positive_above else None,
        'smallest_30': ev[:30].tolist(),
        'zero_proj_frac': float(zero_proj / total_proj)
    }


def main():
    all_results = []

    for N in [32, 48, 64]:
        print(f"\nBuilding N={N}...", flush=True)
        t0 = time.time()
        tau, VM = build_tau(N)
        t1 = time.time()
        print(f"  Built in {t1-t0:.1f}s", flush=True)
        r = analyze(N, tau, VM)
        all_results.append(r)

        with open(OUTFILE, 'w') as f:
            json.dump(all_results, f, indent=2)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()

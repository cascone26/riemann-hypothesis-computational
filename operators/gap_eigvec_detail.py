#!/usr/bin/env python3
"""
gap_eigvec_detail.py — Analyze the eigenvector AT the spectral gap.

The 24th eigenvalue (~1.9e-4) is the first above the near-zero cluster.
This script:
  1. Finds the 24th eigenvector and identifies its dominant n-modes
  2. Computes tau[n,n] for the dominant modes (arch, B, P breakdown)
  3. Analyzes the coupling from n=±12 to near-zero modes
  4. Checks if the gap scales with the distance from the nearest zeta zero
  5. Tracks how the gap eigenvalue changes with N (N=32,48,64,96)
"""

import json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
LAMBDA_SQ = 14
L = float(mp.log(LAMBDA_SQ))
sh = float(np.sinh(L/4)**2)

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/gap_eigvec_detail.json")

# Known zeta zeros
GAMMA = [14.134725, 21.022040, 25.010856, 30.424876, 32.935062]


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
    G = 60000
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
    return (tau + tau.T)/2, arch, B, P, idx


def analyze_gap(N):
    print(f"\n{'='*60}")
    print(f"N={N}: Building tau_N...", flush=True)
    t0 = time.time()
    tau, arch, B, P, idx = build_tau(N)
    print(f"  Built in {time.time()-t0:.1f}s", flush=True)

    ev, evec = np.linalg.eigh(tau)

    # Find near-zero cluster and gap
    nz = int(np.sum(np.abs(ev) < 1e-4))
    above = ev[ev > 1e-4]
    gap_val = float(above[0]) if len(above) > 0 else None

    print(f"  Near-zero count: {nz}")
    print(f"  Spectral gap eigenvalue: {gap_val:.6e}" if gap_val else "  No gap found")

    if gap_val is None:
        return None

    # The gap eigenvector index (first above 1e-4 threshold)
    gap_idx = np.argmax(ev > 1e-4)
    gap_evec = evec[:, gap_idx]  # the gap eigenvector

    # Power per n-mode
    pwr = np.abs(gap_evec)**2
    top8 = np.argsort(pwr)[-8:][::-1]

    print(f"\n  Gap eigenvector (lambda={ev[gap_idx]:.4e}): top modes")
    mode_info = []
    for i in top8:
        n = idx[i]
        s_n = TWO_PI * n / L
        arch_nn = float(arch[i, i])
        B_nn = float(B[i, i])
        P_nn = float(P[i, i])
        tau_nn = float(tau[i, i])
        # Distance from nearest zeta zero
        dists = [(abs(s_n - g), j+1) for j, g in enumerate(GAMMA)]
        d_min, zero_num = min(dists, key=lambda x: x[0])
        print(f"    n={n:>4} ({pwr[i]:.4f}): s_n={s_n:.3f}, tau_nn={tau_nn:.4f} "
              f"(arch={arch_nn:.4f}, B={B_nn:.4f}, P={P_nn:.4f}), "
              f"dist_to_gamma_{zero_num}={d_min:.3f}")
        mode_info.append({'n': n, 'pwr': float(pwr[i]), 's_n': float(s_n),
                          'arch_nn': arch_nn, 'B_nn': B_nn, 'P_nn': P_nn,
                          'tau_nn': tau_nn, 'nearest_zero': int(zero_num),
                          'dist_to_zero': float(d_min)})

    # Also analyze the coupling of n=±12 to lower modes
    i12 = idx.index(12) if 12 in idx else None
    i_m12 = idx.index(-12) if -12 in idx else None

    coupling_info = {}
    if i12 is not None:
        # Row n=12: off-diagonal elements
        row12 = tau[i12, :]
        row12_abs = np.abs(row12)
        row12_abs[i12] = 0.0  # zero out diagonal
        top_coup = np.argsort(row12_abs)[-6:][::-1]
        print(f"\n  Row n=12 coupling (|tau[12,j]| for j≠12):")
        for j in top_coup:
            n_j = idx[j]
            print(f"    to n={n_j:>4}: {row12[j]:.4e}")
        coupling_info = {
            'tau_12_12': float(tau[i12, i12]),
            'gershgorin_radius': float(np.sum(np.abs(tau[i12])) - abs(tau[i12, i12])),
            'top_coupling_modes': [{'n': idx[j], 'coupling': float(row12[j])} for j in top_coup]
        }
        print(f"  tau[12,12] = {tau[i12,i12]:.6f}")
        print(f"  Gershgorin radius R_12 = {coupling_info['gershgorin_radius']:.6f}")
        print(f"  Gershgorin gap = {tau[i12,i12] - coupling_info['gershgorin_radius']:.4e}")

    # Compare with the gap-related mode distance to zeros
    print(f"\n  Distance from s_12={TWO_PI*12/L:.3f} to zeta zeros:")
    for j, g in enumerate(GAMMA):
        print(f"    gamma_{j+1}={g:.3f}: dist={abs(TWO_PI*12/L - g):.3f}")

    return {
        'N': N,
        'n_near_zero': nz,
        'gap_eigenvalue': gap_val,
        'gap_eigvec_top_modes': mode_info,
        'coupling_n12': coupling_info,
        'ev_around_gap': ev[max(0,gap_idx-3):gap_idx+4].tolist()
    }


def main():
    results = []
    for N in [32, 48, 64, 96]:
        r = analyze_gap(N)
        if r:
            results.append(r)
        with open(OUTFILE, 'w') as f:
            json.dump(results, f, indent=2)

    # Summary: how does the gap eigenvalue depend on N?
    print(f"\n{'='*60}")
    print(f"SPECTRAL GAP CONVERGENCE:")
    for r in results:
        print(f"  N={r['N']:>3}: gap={r['gap_eigenvalue']:.6e}, nz={r['n_near_zero']}")

    # Look for convergence
    if len(results) >= 3:
        gaps = [r['gap_eigenvalue'] for r in results]
        Ns = [r['N'] for r in results]
        # fit log-log
        log_N = np.log(Ns)
        log_g = np.log(gaps)
        c = np.polyfit(log_N, log_g, 1)
        print(f"\n  Log-log fit: gap ~ N^{c[0]:.4f} × exp({c[1]:.4f})")
        print(f"  (rate < 0 → gap shrinks with N, rate ≈ 0 → gap converges)")

    print(f"\nResults: {OUTFILE}")


if __name__ == "__main__":
    main()

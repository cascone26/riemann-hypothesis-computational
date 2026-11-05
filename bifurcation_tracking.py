#!/usr/bin/env python3
"""
bifurcation_tracking.py -- Track gap(eps) for tau(eps) = arch + eps*(w02 - prime)
across eps in [0,1], at multiple N values, fine-grained sampling.

Task: Determine if gap(eps) > 0 for all eps, or if there exists eps where gap -> 0 or < 0.
If gap never hits zero, this supports structural simplicity preservation.
Also isolate which component (w02 vs prime) drives the gap collapse.
"""

import json
import time
import numpy as np
from scipy.special import digamma as _digamma
import matplotlib.pyplot as plt
from pathlib import Path

# Constants from connes_decisive_v3
mp_dps = 30
LAMBDA_SQ = 14
L = float(2 * np.log(np.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))

def vm_table(kmax):
    """Von Mangoldt Lambda table for k <= kmax."""
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
            tab.append((k, float(np.log(p))))
    return tab

VM = vm_table(LAMBDA_SQ)

def omega_vec(r):
    """Omega(r) = Re psi(1/4 + i r/2) - log pi, vectorized."""
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI

def e_hat_vec(r, n):
    """e_hat_n(r) vectorized."""
    s = TWO_PI * n / L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L) + 0j
    dd = d[~small]
    out[~small] = (1.0 / np.sqrt(L)) * np.exp(0.5j * dd * L) * 2.0 * np.sin(0.5 * dd * L) / dd
    return out

def q_kernel(n, m, y):
    """q kernel for prime term."""
    if n == m:
        return 2.0 * (1.0 - abs(y) / L) * np.cos(TWO_PI * n * y / L)
    return (np.sin(TWO_PI * m * y / L) - np.sin(TWO_PI * n * y / L)) / (np.pi * (n - m))

def prime_mat(idx):
    """Prime/von Mangoldt matrix."""
    D = len(idx)
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            s = 0.0
            for k, Lam in VM:
                s += Lam * k ** (-0.5) * q_kernel(n, m, np.log(k))
            P[ii, jj] = P[jj, ii] = s
    return P

def w02_mat(idx):
    """W_{0,2} boundary/pole matrix."""
    D = len(idx)
    W = np.zeros((D, D))
    sh = np.sinh(L / 4) ** 2
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = (L ** 2 - 16 * np.pi ** 2 * m * n)
            den = (L ** 2 + 16 * np.pi ** 2 * m ** 2) * (L ** 2 + 16 * np.pi ** 2 * n ** 2)
            W[ii, jj] = 32 * L * sh * num / den
    return W

def arch_mat(N, R, G):
    """Archimedean block."""
    idx = list(range(-N, N + 1))
    D = len(idx)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h)
    w[0] = w[-1] = 0.5 * h
    Om = omega_vec(grid)
    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = e_hat_vec(grid, n)
    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)
    return arch, idx

def build_tau_eps(arch, w02, prime, eps):
    """Construct tau(eps) = arch + eps*(w02 - prime)."""
    pert = w02 - prime
    tau = arch + eps * pert
    tau = (tau + tau.T) / 2  # Ensure symmetry
    return tau

def analyze_spectrum(N, n_eps_steps=100):
    """
    Bifurcation analysis: sweep eps in [0,1] and compute gap(eps).

    Returns:
    - eps_vals: array of eps values
    - gaps: array of spectral gaps (lambda_2 - lambda_1)
    - lambda1_vals: smallest eigenvalue at each eps
    - lambda2_vals: second-smallest eigenvalue at each eps
    - idx: mode indices
    - arch, w02, prime: the component matrices
    """
    print(f"\n{'='*70}")
    print(f"Bifurcation tracking: N={N}, eps steps={n_eps_steps}")
    print(f"{'='*70}")

    # Build component matrices (once, reused across eps)
    R = TWO_PI * N / L + 80.0
    G = max(60000, 3000 * N)
    print(f"N={N}, dim={2*N+1}, R={R:.1f}, grid_size={G}")

    t0 = time.time()
    arch, idx = arch_mat(N, R, G)
    print(f"  Computed arch: {time.time()-t0:.2f}s")

    t0 = time.time()
    w02 = w02_mat(idx)
    print(f"  Computed w02: {time.time()-t0:.2f}s")

    t0 = time.time()
    prime = prime_mat(idx)
    print(f"  Computed prime: {time.time()-t0:.2f}s")

    # Compute norms for reference
    norm_arch = np.linalg.norm(arch, 'fro')
    norm_pert = np.linalg.norm(w02 - prime, 'fro')
    print(f"  ‖arch‖_F = {norm_arch:.4f}")
    print(f"  ‖w02-prime‖_F = {norm_pert:.4f}")
    print(f"  Ratio = {norm_pert/norm_arch:.4f}")

    # Sweep eps
    eps_vals = np.linspace(0, 1, n_eps_steps)
    gaps = []
    lambda1_vals = []
    lambda2_vals = []

    print(f"\nSweeping eps from 0 to 1 ({n_eps_steps} steps)...")
    t_start = time.time()
    for i, eps in enumerate(eps_vals):
        tau = build_tau_eps(arch, w02, prime, eps)
        evals = np.linalg.eigvalsh(tau)  # Sorted ascending
        gap = evals[1] - evals[0]
        gaps.append(gap)
        lambda1_vals.append(evals[0])
        lambda2_vals.append(evals[1])

        if i % max(1, n_eps_steps // 10) == 0:
            print(f"  eps={eps:.3f}: λ₁={evals[0]:+.4e}, λ₂={evals[1]:+.4e}, gap={gap:+.4e}")

    elapsed = time.time() - t_start
    print(f"Total sweep time: {elapsed:.2f}s")

    gaps = np.array(gaps)
    lambda1_vals = np.array(lambda1_vals)
    lambda2_vals = np.array(lambda2_vals)

    # Analysis
    min_gap = np.min(gaps)
    max_gap = np.max(gaps)
    min_gap_idx = np.argmin(gaps)
    eps_min_gap = eps_vals[min_gap_idx]

    print(f"\n--- Gap Analysis ---")
    print(f"Max gap: {max_gap:.4e} at eps={eps_vals[np.argmax(gaps)]:.3f}")
    print(f"Min gap: {min_gap:.4e} at eps={eps_min_gap:.3f}")
    print(f"Gap ratio (max/min): {max_gap/abs(min_gap) if min_gap != 0 else np.inf:.2e}")

    # Check if gap crosses zero
    sign_changes = np.sum(np.diff(np.sign(gaps)) != 0)
    print(f"Zero crossings: {sign_changes}")

    if min_gap < 0:
        print(f"⚠️  WARNING: Gap went NEGATIVE! Min gap = {min_gap:.4e}")
        crossings = np.where(gaps < 0)[0]
        print(f"   Negative at {len(crossings)} points")
        for j in crossings[:10]:  # Show first 10
            print(f"     eps={eps_vals[j]:.4f}: gap={gaps[j]:.4e}")
    else:
        print(f"✓ Gap stays POSITIVE throughout. Min gap = {min_gap:.4e}")

    return {
        'N': N,
        'eps_vals': eps_vals,
        'gaps': gaps,
        'lambda1_vals': lambda1_vals,
        'lambda2_vals': lambda2_vals,
        'idx': idx,
        'arch': arch,
        'w02': w02,
        'prime': prime,
        'min_gap': min_gap,
        'max_gap': max_gap,
        'eps_min_gap': eps_min_gap,
    }

def analyze_components(N, n_eps_steps=50):
    """
    Isolate the perturbation components: analyze tau_w02 (arch+w02) and tau_prime (arch-prime) separately.
    """
    print(f"\n{'='*70}")
    print(f"Component Isolation: N={N}")
    print(f"{'='*70}")

    # Build components
    R = TWO_PI * N / L + 80.0
    G = max(60000, 3000 * N)
    arch, idx = arch_mat(N, R, G)
    w02 = w02_mat(idx)
    prime = prime_mat(idx)

    # Test 1: arch + w02 (no prime correction)
    eps_vals = np.linspace(0, 1, n_eps_steps)
    gaps_w02 = []

    print(f"\nTest 1: tau_w02(eps) = arch + eps*w02 (no prime)")
    for eps in eps_vals:
        tau = arch + eps * w02
        tau = (tau + tau.T) / 2
        evals = np.linalg.eigvalsh(tau)
        gap = evals[1] - evals[0]
        gaps_w02.append(gap)

    min_gap_w02 = np.min(gaps_w02)
    print(f"  Min gap: {min_gap_w02:.4e}")
    print(f"  Gap stays positive? {min_gap_w02 > 0}")

    # Test 2: arch - prime (no w02 correction)
    gaps_prime = []

    print(f"\nTest 2: tau_prime(eps) = arch - eps*prime (no w02)")
    for eps in eps_vals:
        tau = arch - eps * prime
        tau = (tau + tau.T) / 2
        evals = np.linalg.eigvalsh(tau)
        gap = evals[1] - evals[0]
        gaps_prime.append(gap)

    min_gap_prime = np.min(gaps_prime)
    print(f"  Min gap: {min_gap_prime:.4e}")
    print(f"  Gap stays positive? {min_gap_prime > 0}")

    # Test 3: Full bifurcation (w02 and prime together)
    gaps_full = []

    print(f"\nTest 3: tau(eps) = arch + eps*(w02-prime) (full)")
    for eps in eps_vals:
        pert = w02 - prime
        tau = arch + eps * pert
        tau = (tau + tau.T) / 2
        evals = np.linalg.eigvalsh(tau)
        gap = evals[1] - evals[0]
        gaps_full.append(gap)

    min_gap_full = np.min(gaps_full)
    print(f"  Min gap: {min_gap_full:.4e}")
    print(f"  Gap stays positive? {min_gap_full > 0}")

    return {
        'N': N,
        'eps_vals': eps_vals,
        'gaps_w02': np.array(gaps_w02),
        'gaps_prime': np.array(gaps_prime),
        'gaps_full': np.array(gaps_full),
        'min_gap_w02': min_gap_w02,
        'min_gap_prime': min_gap_prime,
        'min_gap_full': min_gap_full,
    }

def plot_results(results_list, comp_results_list):
    """Plot gap vs eps for multiple N."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Gap vs eps for different N (full bifurcation)
    ax = axes[0, 0]
    for res in results_list:
        ax.semilogy(res['eps_vals'], np.abs(res['gaps']), label=f"N={res['N']}", marker='o', markersize=3)
    ax.axhline(y=0, color='r', linestyle='--', label='zero')
    ax.set_xlabel('ε')
    ax.set_ylabel('|Gap(ε)| = |λ₂ - λ₁|')
    ax.set_title('Bifurcation Tracking: Full tau(ε) = arch + ε(w02-prime)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Component isolation (single N)
    if comp_results_list:
        ax = axes[0, 1]
        comp = comp_results_list[0]
        ax.semilogy(comp['eps_vals'], np.abs(comp['gaps_w02']), label='w02 only', marker='s', markersize=3)
        ax.semilogy(comp['eps_vals'], np.abs(comp['gaps_prime']), label='prime only', marker='^', markersize=3)
        ax.semilogy(comp['eps_vals'], np.abs(comp['gaps_full']), label='w02-prime', marker='o', markersize=3)
        ax.axhline(y=0, color='r', linestyle='--', alpha=0.5)
        ax.set_xlabel('ε')
        ax.set_ylabel('|Gap(ε)|')
        ax.set_title(f'Component Isolation: N={comp["N"]}')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Plot 3: λ₁ vs eps
    ax = axes[1, 0]
    for res in results_list:
        ax.plot(res['eps_vals'], res['lambda1_vals'], label=f"N={res['N']}", marker='o', markersize=3)
    ax.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax.set_xlabel('ε')
    ax.set_ylabel('λ₁')
    ax.set_title('Smallest Eigenvalue vs ε')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: λ₂ vs eps
    ax = axes[1, 1]
    for res in results_list:
        ax.plot(res['eps_vals'], res['lambda2_vals'], label=f"N={res['N']}", marker='o', markersize=3)
    ax.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax.set_xlabel('ε')
    ax.set_ylabel('λ₂')
    ax.set_title('Second-Smallest Eigenvalue vs ε')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/Users/scones/projects/riemann/results/bifurcation_tracking.png', dpi=150)
    print(f"\nPlot saved to results/bifurcation_tracking.png")
    plt.close()

def main():
    print("\n" + "="*70)
    print("BIFURCATION TRACKING: tau(eps) = arch + eps(w02 - prime)")
    print("="*70)

    # Run at multiple N values with fine-grained eps sampling
    N_values = [24, 32]  # Two representative sizes
    results_list = []

    for N in N_values:
        res = analyze_spectrum(N, n_eps_steps=100)
        results_list.append(res)

    # Component isolation at one N
    comp_results_list = []
    comp_res = analyze_components(N=32, n_eps_steps=50)
    comp_results_list.append(comp_res)

    # Plot
    plot_results(results_list, comp_results_list)

    # Save results
    output_file = Path('/Users/scones/projects/riemann/results/bifurcation_tracking.json')
    output_file.parent.mkdir(parents=True, exist_ok=True)

    save_data = {
        'timestamp': str(time.time()),
        'summary': {
            'task': 'Bifurcation tracking: tau(eps) = arch + eps*(w02-prime), eps in [0,1]',
            'n_eps_steps': 100,
            'N_values': N_values,
        },
        'results': [
            {
                'N': res['N'],
                'eps_vals': res['eps_vals'].tolist(),
                'gaps': res['gaps'].tolist(),
                'lambda1_vals': res['lambda1_vals'].tolist(),
                'lambda2_vals': res['lambda2_vals'].tolist(),
                'min_gap': float(res['min_gap']),
                'max_gap': float(res['max_gap']),
                'eps_min_gap': float(res['eps_min_gap']),
            }
            for res in results_list
        ],
        'component_analysis': [
            {
                'N': comp['N'],
                'eps_vals': comp['eps_vals'].tolist(),
                'gaps_w02': comp['gaps_w02'].tolist(),
                'gaps_prime': comp['gaps_prime'].tolist(),
                'gaps_full': comp['gaps_full'].tolist(),
                'min_gap_w02': float(comp['min_gap_w02']),
                'min_gap_prime': float(comp['min_gap_prime']),
                'min_gap_full': float(comp['min_gap_full']),
            }
            for comp in comp_results_list
        ]
    }

    with open(output_file, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\n✓ Results saved to {output_file}")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    for res in results_list:
        print(f"\nN={res['N']}:")
        print(f"  Min gap (across eps=[0,1]): {res['min_gap']:.4e}")
        print(f"  Max gap: {res['max_gap']:.4e}")
        print(f"  Min at eps: {res['eps_min_gap']:.3f}")
        print(f"  Gap stays positive? {'YES ✓' if res['min_gap'] > 0 else 'NO ✗'}")

    print(f"\n{'='*70}")
    print("COMPONENT ANALYSIS (N=32)")
    print(f"{'='*70}")
    if comp_results_list:
        comp = comp_results_list[0]
        print(f"w02 only (arch + eps*w02):")
        print(f"  Min gap: {comp['min_gap_w02']:.4e}, stays positive? {comp['min_gap_w02'] > 0}")
        print(f"prime only (arch - eps*prime):")
        print(f"  Min gap: {comp['min_gap_prime']:.4e}, stays positive? {comp['min_gap_prime'] > 0}")
        print(f"Full (arch + eps*(w02-prime)):")
        print(f"  Min gap: {comp['min_gap_full']:.4e}, stays positive? {comp['min_gap_full'] > 0}")

if __name__ == '__main__':
    main()

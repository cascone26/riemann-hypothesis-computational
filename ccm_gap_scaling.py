#!/usr/bin/env python3
"""
ccm_gap_scaling.py -- Compute minimum eigenvalue gap for tau(eps=1) across N values.

Task: Extend bifurcation_tracking findings (N=24,32 with 100-point eps sweep) to a
WIDER range of N (8,12,16,20,24,28,32,40,48 and more if tractable), computing ONLY
the full-operator (eps=1) minimum gap at each N, and then fitting scaling laws.

Returns:
- gap(N) data table for each computed N
- Fitted scaling law (power law, exponential, or other) with parameters + R²
- Sanity check: do no sign changes appear at larger N?
"""

import json
import time
import numpy as np
from scipy.special import digamma as _digamma
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pathlib import Path

# Constants from connes_decisive_v3
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
    tau = (tau + tau.T) / 2
    return tau

def compute_gap_at_N(N, eps=1.0, n_eps_check=5):
    """
    Compute the gap at a specific N and eps.

    Also do a cheap eps-sweep sanity check: evaluate gap at n_eps_check points
    in [0, eps] to check for sign changes or anomalies.

    Returns: (gap_at_eps, min_gap_in_range, max_gap_in_range, eps_check_vals, gap_check_vals, success)
    """
    print(f"\n{'='*70}")
    print(f"Computing gap for N={N}")
    print(f"{'='*70}")

    # Build components
    R = TWO_PI * N / L + 80.0
    G = max(60000, 3000 * N)
    print(f"N={N}, dim={2*N+1}, R={R:.1f}, grid_size={G}")

    t0 = time.time()
    try:
        arch, idx = arch_mat(N, R, G)
        print(f"  Computed arch: {time.time()-t0:.2f}s")

        t0 = time.time()
        w02 = w02_mat(idx)
        print(f"  Computed w02: {time.time()-t0:.2f}s")

        t0 = time.time()
        prime = prime_mat(idx)
        print(f"  Computed prime: {time.time()-t0:.2f}s")
    except Exception as e:
        print(f"  ERROR building components: {e}")
        return None, None, None, [], [], False

    # Compute gap at full eps
    t0 = time.time()
    try:
        tau = build_tau_eps(arch, w02, prime, eps)
        evals = np.linalg.eigvalsh(tau)
        gap = evals[1] - evals[0]
        gap_time = time.time() - t0
        print(f"  Gap at eps={eps}: {gap:.6e} ({gap_time:.2f}s)")
        print(f"  λ₁={evals[0]:.6e}, λ₂={evals[1]:.6e}")
    except Exception as e:
        print(f"  ERROR computing eigenvalues at eps={eps}: {e}")
        return None, None, None, [], [], False

    # Sanity check: sample a few eps values and look for anomalies
    eps_check_vals = np.linspace(0, eps, n_eps_check)
    gap_check_vals = []

    print(f"\n  Sanity check: sampling gap at {n_eps_check} eps points...")
    for eps_i in eps_check_vals:
        try:
            tau_i = build_tau_eps(arch, w02, prime, eps_i)
            evals_i = np.linalg.eigvalsh(tau_i)
            gap_i = evals_i[1] - evals_i[0]
            gap_check_vals.append(gap_i)
            status = "✓" if gap_i > 0 else "✗ NEGATIVE"
            print(f"    eps={eps_i:.3f}: gap={gap_i:.6e} {status}")
        except Exception as e:
            print(f"    eps={eps_i:.3f}: ERROR: {e}")
            gap_check_vals.append(None)

    gap_check_vals = np.array(gap_check_vals)

    # Check for sign changes
    valid_gaps = gap_check_vals[gap_check_vals != None]
    if len(valid_gaps) > 0:
        min_gap_in_range = float(np.min(valid_gaps))
        max_gap_in_range = float(np.max(valid_gaps))
        if np.any(valid_gaps < 0):
            print(f"  ⚠️  WARNING: gap went NEGATIVE in check!")
        else:
            print(f"  ✓ Gap stays POSITIVE in check")
    else:
        min_gap_in_range = None
        max_gap_in_range = None

    return float(gap), min_gap_in_range, max_gap_in_range, eps_check_vals.tolist(), gap_check_vals.tolist(), True

def fit_scaling_laws(N_vals, gaps):
    """
    Fit gap(N) to candidate scaling laws:
    - Power law: gap(N) = C * N^(-alpha)
    - Exponential: gap(N) = C * exp(-beta * N)
    - Mixed: gap(N) = C * N^(-alpha) * exp(-beta * N)

    Returns dict with fit results, parameters, and R².
    """
    N_vals = np.array(N_vals, dtype=float)
    gaps = np.array(gaps, dtype=float)

    results = {}

    # Filter out invalid gaps
    valid = gaps > 0
    if not np.any(valid):
        print("No valid gap data to fit.")
        return results

    N_valid = N_vals[valid]
    gaps_valid = gaps[valid]

    print(f"\n{'='*70}")
    print("FITTING SCALING LAWS")
    print(f"{'='*70}")
    print(f"Valid data points: {len(N_valid)}")
    print(f"N range: [{N_valid.min():.0f}, {N_valid.max():.0f}]")
    print(f"Gap range: [{gaps_valid.min():.2e}, {gaps_valid.max():.2e}]")

    # Power law: log(gap) = log(C) - alpha*log(N)
    try:
        def power_law(N, C, alpha):
            return C * N ** (-alpha)

        popt, _ = curve_fit(power_law, N_valid, gaps_valid, p0=[1e-6, 1.0], maxfev=10000)
        C_power, alpha_power = popt
        gaps_fit = power_law(N_valid, C_power, alpha_power)
        ss_res = np.sum((gaps_valid - gaps_fit) ** 2)
        ss_tot = np.sum((gaps_valid - np.mean(gaps_valid)) ** 2)
        r2_power = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        results['power_law'] = {
            'formula': 'gap(N) = C * N^(-alpha)',
            'C': float(C_power),
            'alpha': float(alpha_power),
            'R2': float(r2_power),
        }
        print(f"\nPower Law: gap(N) = {C_power:.4e} * N^(-{alpha_power:.4f})")
        print(f"  R² = {r2_power:.6f}")
    except Exception as e:
        print(f"\nPower law fit failed: {e}")

    # Exponential: log(gap) = log(C) - beta*N
    try:
        def exp_law(N, C, beta):
            return C * np.exp(-beta * N)

        popt, _ = curve_fit(exp_law, N_valid, gaps_valid, p0=[1e-1, 0.1], maxfev=10000)
        C_exp, beta_exp = popt
        gaps_fit = exp_law(N_valid, C_exp, beta_exp)
        ss_res = np.sum((gaps_valid - gaps_fit) ** 2)
        ss_tot = np.sum((gaps_valid - np.mean(gaps_valid)) ** 2)
        r2_exp = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        results['exponential'] = {
            'formula': 'gap(N) = C * exp(-beta*N)',
            'C': float(C_exp),
            'beta': float(beta_exp),
            'R2': float(r2_exp),
        }
        print(f"\nExponential: gap(N) = {C_exp:.4e} * exp(-{beta_exp:.6f}*N)")
        print(f"  R² = {r2_exp:.6f}")
    except Exception as e:
        print(f"\nExponential fit failed: {e}")

    # Mixed: log(gap) = log(C) - alpha*log(N) - beta*N
    try:
        def mixed_law(N, C, alpha, beta):
            return C * N ** (-alpha) * np.exp(-beta * N)

        popt, _ = curve_fit(mixed_law, N_valid, gaps_valid, p0=[1e-3, 0.5, 0.01], maxfev=10000)
        C_mix, alpha_mix, beta_mix = popt
        gaps_fit = mixed_law(N_valid, C_mix, alpha_mix, beta_mix)
        ss_res = np.sum((gaps_valid - gaps_fit) ** 2)
        ss_tot = np.sum((gaps_valid - np.mean(gaps_valid)) ** 2)
        r2_mix = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        results['mixed'] = {
            'formula': 'gap(N) = C * N^(-alpha) * exp(-beta*N)',
            'C': float(C_mix),
            'alpha': float(alpha_mix),
            'beta': float(beta_mix),
            'R2': float(r2_mix),
        }
        print(f"\nMixed: gap(N) = {C_mix:.4e} * N^(-{alpha_mix:.4f}) * exp(-{beta_mix:.6f}*N)")
        print(f"  R² = {r2_mix:.6f}")
    except Exception as e:
        print(f"\nMixed fit failed: {e}")

    return results

def plot_results(N_vals, gaps, fits):
    """Plot gap(N) data and fitted curves."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Data points
    valid = np.array(gaps) > 0
    N_valid = np.array(N_vals)[valid]
    gaps_valid = np.array(gaps)[valid]
    ax.semilogy(N_valid, gaps_valid, 'ko', markersize=8, label='Computed gap')

    # Fitted curves
    N_smooth = np.linspace(N_valid.min(), N_valid.max(), 200)

    if 'power_law' in fits:
        fit = fits['power_law']
        C, alpha = fit['C'], fit['alpha']
        gaps_fit = C * N_smooth ** (-alpha)
        ax.semilogy(N_smooth, gaps_fit, 'b-', linewidth=2,
                   label=f"Power: C*N^(-α), α={alpha:.3f}, R²={fit['R2']:.4f}")

    if 'exponential' in fits:
        fit = fits['exponential']
        C, beta = fit['C'], fit['beta']
        gaps_fit = C * np.exp(-beta * N_smooth)
        ax.semilogy(N_smooth, gaps_fit, 'r-', linewidth=2,
                   label=f"Exp: C*exp(-β*N), β={beta:.6f}, R²={fit['R2']:.4f}")

    if 'mixed' in fits:
        fit = fits['mixed']
        C, alpha, beta = fit['C'], fit['alpha'], fit['beta']
        gaps_fit = C * N_smooth ** (-alpha) * np.exp(-beta * N_smooth)
        ax.semilogy(N_smooth, gaps_fit, 'g--', linewidth=2,
                   label=f"Mixed: C*N^(-α)*exp(-β*N), α={alpha:.3f}, β={beta:.6f}, R²={fit['R2']:.4f}")

    ax.set_xlabel('N (bandwidth)', fontsize=12)
    ax.set_ylabel('min gap(λ₂ - λ₁) at eps=1', fontsize=12)
    ax.set_title('CCM Operator: Eigenvalue Gap Scaling', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    outfile = '/Users/scones/projects/riemann/results/ccm_gap_scaling.png'
    plt.savefig(outfile, dpi=150)
    print(f"\nPlot saved to {outfile}")
    plt.close()

def main():
    print("\n" + "="*70)
    print("CCM OPERATOR: MINIMUM EIGENVALUE GAP SCALING ANALYSIS")
    print("="*70)
    print("Task: compute gap(N) for tau(eps=1) across a range of N,")
    print("      then fit scaling laws and check for anomalies at large N.")

    # N values to compute: aim for coverage, skip expensive ones if they hit timeout
    N_target = [8, 12, 16, 20, 24, 28, 32, 40, 48]
    N_vals = []
    gaps = []
    metadata = []

    max_time_per_N = 300  # 5 min per N as timeout

    for N in N_target:
        start = time.time()
        gap, min_gap_check, max_gap_check, eps_check, gap_check, success = compute_gap_at_N(N, eps=1.0, n_eps_check=5)
        elapsed = time.time() - start

        if success and gap is not None:
            N_vals.append(N)
            gaps.append(gap)
            metadata.append({
                'N': N,
                'gap': gap,
                'min_gap_in_check': min_gap_check,
                'max_gap_in_check': max_gap_check,
                'eps_check': eps_check,
                'gap_check': gap_check,
                'time_sec': elapsed,
            })
            print(f"\n✓ N={N}: gap={gap:.6e}, time={elapsed:.1f}s")

            # If N is getting slow, stop early
            if elapsed > max_time_per_N:
                print(f"  ⚠️  Taking >5min per N; stopping here for tractability.")
                break
        else:
            print(f"\n✗ N={N}: FAILED (or timeout)")

    # Fit scaling laws
    fits = fit_scaling_laws(N_vals, gaps)

    # Plot
    if len(N_vals) > 1:
        plot_results(N_vals, gaps, fits)

    # Save results
    output_file = Path('/Users/scones/projects/riemann/results/ccm_gap_scaling.json')
    output_file.parent.mkdir(parents=True, exist_ok=True)

    save_data = {
        'timestamp': str(time.time()),
        'summary': {
            'task': 'CCM operator gap(N) scaling: compute min eigenvalue gap at eps=1 across N',
            'N_target': N_target,
            'N_computed': N_vals,
        },
        'data': metadata,
        'scaling_fits': fits,
    }

    with open(output_file, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\n✓ Results saved to {output_file}")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"\nComputed {len(N_vals)} N values:")
    for n, g in zip(N_vals, gaps):
        print(f"  N={n:2d}: gap = {g:.6e}")

    print(f"\nScaling law fits:")
    for name, fit in fits.items():
        print(f"\n  {name}:")
        print(f"    Formula: {fit['formula']}")
        for key, val in fit.items():
            if key != 'formula':
                if isinstance(val, float):
                    print(f"    {key} = {val:.6e}" if abs(val) < 1e-3 else f"    {key} = {val:.6f}")
                else:
                    print(f"    {key} = {val}")

if __name__ == '__main__':
    main()

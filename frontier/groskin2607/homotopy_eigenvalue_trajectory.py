#!/usr/bin/env python3
"""
Eigenvalue trajectory analysis via homotopy continuation.

Track how λ_min(t) evolves as we deform τ(t) = (W02-WR) - t·Wp from t=0 (background only)
to t=1 (full matrix with primes). Look for avoided crossings and level repulsion that might
explain why the final margin is so small and scales like c^(-0.44N).

Key idea: The full coupling (S2642) means eigenvalues interact globally. Even if individual
perturbations are small, near-degeneracies at some intermediate t can cause eigenvalues to
repel sharply, pushing one down toward zero.

Usage:
  python3 homotopy_eigenvalue_trajectory.py --c 30 --N 8 --steps 51 --plot

Generates:
  - trajectory data (eigenvalue vs t for several small/medium eigenvalues)
  - plots of level crossings and avoided crossings
  - estimates of crossing point, repulsion distance, local scaling
"""

import sys
import argparse
import numpy as np
import json
from pathlib import Path

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from arb_ldlt_certify import build_arb_tau

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from flint import arb, ctx
import mpmath as mp


def build_tau_float64(c, N, t_weight):
    """
    Build τ(t) = (W02 - WR) - t_weight * Wp at float64 precision.
    Returns the matrix and the dimension.

    Args:
        c: cutoff parameter
        N: number of frequency modes
        t_weight: parameter t ∈ [0,1] controlling how much primes to include
                  (t=0: background only, t=1: full matrix)
    """
    dim = 2*N + 1
    tau = np.zeros((dim, dim))

    # Build components separately using mpmath for intermediate precision
    mp.mp.dps = 50  # 50 decimal digits ≈ 166 bits
    from arb_ldlt_certify import prime_powers_up_to, arb_closed_forms, arb_J, arb_kappa

    # Use Arb for precision, then convert to float64
    ctx.prec = 300
    S_arb, CC_arb, XC_arb, L_arb = arb_closed_forms(N, c, 300)

    # Convert to mpmath
    L = mp.log(c)
    S = [float(s) if s != 0 else 0.0 for s in S_arb]
    CC = [float(cc) if cc != 0 else 0.0 for cc in CC_arb]
    XC = [float(xc) if xc != 0 else 0.0 for xc in XC_arb]

    PI = mp.pi
    sp2 = 16 * PI * PI
    l2 = L * L
    pref02 = 32 * L * mp.sinh(L/4) ** 2

    # Get kappa and J
    eL = mp.exp(L)
    kappa = mp.log(4 * PI * (eL - 1) / (eL + 1)) + mp.euler
    U = mp.exp(L/2)
    J = -2*mp.log(U+1) + mp.log(U*U + 1) + 2*mp.atan(U) + mp.log(2) - PI/2

    # Prime data
    pdata = prime_powers_up_to(c)
    weights = [mp.log(p) / mp.sqrt(q) for (q, p) in pdata]
    positions = [mp.log(q) for (q, p) in pdata]
    M = len(weights)

    def S_signed(nn):
        return S[nn] if nn >= 0 else -S[-nn]

    # Build matrix
    for i in range(dim):
        n = i - N
        for j in range(i, dim):
            m = j - N

            # W02 term
            num = l2 - sp2 * m * n
            den = (l2 + sp2*m*m) * (l2 + sp2*n*n)
            W02 = float(pref02 * num / den)

            # WR term
            if n == m:
                WR = float(kappa + 2*CC[abs(n)] + J - (2/L)*XC[abs(n)])
            else:
                WR = float((S_signed(m) - S_signed(n)) / (PI * (n - m)))

            # Wp term (primes)
            Wp = 0.0
            for idx in range(M):
                y = positions[idx]
                if n == m:
                    q = 2 * (1 - y/L) * mp.cos(2*PI*n*y/L)
                else:
                    q = (mp.sin(2*PI*m*y/L) - mp.sin(2*PI*n*y/L)) / (PI*(n-m))
                Wp += float(weights[idx] * q)

            # τ(t) = (W02 - WR) - t * Wp
            val = W02 - WR - t_weight * Wp
            tau[i,j] = val
            if i != j:
                tau[j,i] = val

    return tau, dim


def track_eigenvalue_trajectory(c, N, num_steps=51, focus_k=None):
    """
    Compute the trajectory of the smallest eigenvalues as t goes from 0 to 1.

    Args:
        c, N: matrix parameters
        num_steps: number of t values (default 51 gives t=0, 0.02, 0.04, ..., 1.0)
        focus_k: if specified, track the k-th smallest eigenvalues (not just the global min)
                 use focus_k=5 to track the 5 smallest

    Returns:
        dict with:
          - t_values: array of t parameters
          - eigenvalue_trajectories: dict mapping eigenvalue_index -> (t_array, lambda_array)
          - metadata: (c, N, num_steps)
    """
    if focus_k is None:
        focus_k = min(10, 2*N+1)  # track up to 10 eigenvalues or all if dim is small

    t_values = np.linspace(0, 1, num_steps)
    trajectories = {}

    print(f"Tracking eigenvalue trajectories for c={c}, N={N}:")
    print(f"  dim = {2*N+1}, focusing on {focus_k} smallest eigenvalues")
    print(f"  t ∈ [0,1], {num_steps} steps")

    for step_idx, t in enumerate(t_values):
        tau, dim = build_tau_float64(c, N, t)

        # Compute eigenvalues (sorted ascending)
        try:
            evals = np.linalg.eigvalsh(tau)
            evals_sorted = np.sort(evals)
        except np.linalg.LinAlgError:
            print(f"  WARNING: eigenvalue computation failed at t={t:.3f}, skipping")
            continue

        # Store the k smallest
        for k in range(min(focus_k, len(evals_sorted))):
            if k not in trajectories:
                trajectories[k] = ([], [])
            trajectories[k][0].append(t)
            trajectories[k][1].append(evals_sorted[k])

        if (step_idx + 1) % max(1, num_steps // 5) == 0 or step_idx == 0:
            print(f"  step {step_idx+1:2d}/{num_steps}: t={t:.2f}", end="")
            if len(evals_sorted) >= 1:
                print(f"  λ₀={evals_sorted[0]:+.4e}  λ₁={evals_sorted[1]:+.4e}", end="")
            print()

    # Convert to numpy arrays for easier analysis
    for k in trajectories:
        trajectories[k] = (np.array(trajectories[k][0]), np.array(trajectories[k][1]))

    return {
        't_values': t_values,
        'trajectories': trajectories,
        'metadata': {'c': c, 'N': N, 'num_steps': num_steps},
    }


def analyze_avoided_crossings(result, tol=1e-3):
    """
    Look for avoided crossings in the trajectory data.

    An avoided crossing is detected when:
    1. Two eigenvalue trajectories get close (gap < tol)
    2. Their order swaps (or nearly swaps)
    3. Local curvature is high (indicating repulsion)
    """
    trajectories = result['trajectories']
    c, N, num_steps = result['metadata']['c'], result['metadata']['N'], result['metadata']['num_steps']

    print(f"\n" + "="*100)
    print(f"AVOIDED CROSSING ANALYSIS (c={c}, N={N})")
    print("="*100)

    crossings = []

    # Look for pairs of nearby trajectories
    keys = sorted(trajectories.keys())
    for idx_k in range(len(keys) - 1):
        k1, k2 = keys[idx_k], keys[idx_k + 1]
        t_arr1, lam1 = trajectories[k1]
        t_arr2, lam2 = trajectories[k2]

        if len(t_arr1) < 3 or len(t_arr2) < 3:
            continue

        # Find minimum gap
        min_gap_idx = np.argmin(np.abs(lam1 - lam2))
        min_gap = np.abs(lam1[min_gap_idx] - lam2[min_gap_idx])
        t_gap = t_arr1[min_gap_idx]

        if min_gap < tol and min_gap > 0:
            # Compute local "curvature" = second derivative
            if min_gap_idx > 0 and min_gap_idx < len(t_arr1) - 1:
                dt = t_arr1[min_gap_idx+1] - t_arr1[min_gap_idx-1]
                d2_lam1 = (lam1[min_gap_idx+1] - 2*lam1[min_gap_idx] + lam1[min_gap_idx-1]) / (dt*dt)
                d2_lam2 = (lam2[min_gap_idx+1] - 2*lam2[min_gap_idx] + lam2[min_gap_idx-1]) / (dt*dt)
                curvature = abs(d2_lam1 - d2_lam2)
            else:
                curvature = None

            crossings.append({
                'k1': k1, 'k2': k2,
                't_cross': t_gap,
                'min_gap': min_gap,
                'curvature': curvature,
                'lam1_at_cross': lam1[min_gap_idx],
                'lam2_at_cross': lam2[min_gap_idx],
            })

    # Sort by minimum gap
    crossings.sort(key=lambda x: x['min_gap'])

    for crossing in crossings[:10]:  # Show top 10
        print(f"\nCrossing between λ_{crossing['k1']} and λ_{crossing['k2']}:")
        print(f"  t_cross = {crossing['t_cross']:.4f}")
        print(f"  min_gap = {crossing['min_gap']:.4e}")
        print(f"  λ values at crossing: {crossing['lam1_at_cross']:+.4e}, {crossing['lam2_at_cross']:+.4e}")
        if crossing['curvature'] is not None:
            print(f"  local curvature = {crossing['curvature']:.4e}")

    return crossings


def extract_scaling_laws_from_trajectory(result):
    """
    From the trajectory data, extract local scaling laws:
    - How does λ₀(t) scale near t=1?
    - Is there a power-law or exponential scaling?
    """
    trajectories = result['trajectories']
    c, N = result['metadata']['c'], result['metadata']['N']

    print(f"\n" + "="*100)
    print(f"SCALING LAWS FROM TRAJECTORY (c={c}, N={N})")
    print("="*100)

    if 0 in trajectories:
        t_arr, lam = trajectories[0]

        # Focus on the final portion (t > 0.5) where the action happens
        mask = t_arr > 0.5
        t_final = t_arr[mask]
        lam_final = lam[mask]

        # Avoid taking log of negative/zero values
        lam_final_abs = np.abs(lam_final)
        lam_final_abs[lam_final_abs < 1e-15] = 1e-15

        print(f"\nSmallest eigenvalue λ₀(t) from t=0.5 to t=1:")
        for t_val, lam_val in zip(t_final, lam_final):
            print(f"  t={t_val:.2f}  λ₀={lam_val:+.6e}  log₁₀|λ₀|={np.log10(abs(lam_val) + 1e-20):+.2f}")

        # Try power-law fit in final regime
        if len(t_final) >= 3 and np.all(lam_final > 0):
            log_t = np.log(1 - t_final)  # closer to t=1 means larger (1-t)
            log_lam = np.log(lam_final)

            # Linear fit: log(λ) = a*log(1-t) + b
            coeffs = np.polyfit(log_t, log_lam, 1)
            a_pow, b_pow = coeffs[0], coeffs[1]
            print(f"\nPower-law fit in final regime (t > 0.5):")
            print(f"  log(λ) = {a_pow:.3f}*log(1-t) + {b_pow:.3f}")
            print(f"  → λ ~ (1-t)^{a_pow:.3f}")

            # Check residual
            log_lam_fit = np.polyval(coeffs, log_t)
            resid = np.std(log_lam - log_lam_fit)
            print(f"  fit residual (stdev): {resid:.4f}")


def plot_trajectories(result, output_file=None):
    """
    Plot the eigenvalue trajectories.
    """
    if not HAS_MATPLOTLIB:
        print("(matplotlib not available, skipping plot)")
        return

    trajectories = result['trajectories']
    c, N = result['metadata']['c'], result['metadata']['N']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Left: log scale
    keys = sorted(trajectories.keys())[:8]  # up to 8 eigenvalues
    colors = plt.cm.viridis(np.linspace(0, 1, len(keys)))

    for idx, k in enumerate(keys):
        t_arr, lam = trajectories[k]
        lam_abs = np.abs(lam)
        lam_abs[lam_abs < 1e-15] = 1e-15

        ax1.semilogy(t_arr, lam_abs, 'o-', color=colors[idx], label=f'λ_{k}', markersize=4)

    ax1.set_xlabel('t (weight of primes)', fontsize=12)
    ax1.set_ylabel('|λ| (log scale)', fontsize=12)
    ax1.set_title(f'Eigenvalue trajectories (c={c}, N={N})', fontsize=13)
    ax1.legend(fontsize=9, loc='best')
    ax1.grid(True, alpha=0.3)

    # Right: linear scale, focusing on smallest eigenvalues
    for idx, k in enumerate(keys[:5]):
        t_arr, lam = trajectories[k]
        ax2.plot(t_arr, lam, 'o-', color=colors[idx], label=f'λ_{k}', markersize=4)

    ax2.set_xlabel('t (weight of primes)', fontsize=12)
    ax2.set_ylabel('λ (linear scale)', fontsize=12)
    ax2.set_title(f'Smallest eigenvalues zoomed (c={c}, N={N})', fontsize=13)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(0, color='k', linestyle='--', alpha=0.2)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=120)
        print(f"Saved plot to {output_file}")
    else:
        plt.savefig(f'homotopy_trajectory_c{c}_N{N}.png')
        print(f"Saved plot to homotopy_trajectory_c{c}_N{N}.png")
    plt.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--c', type=int, default=30, help='cutoff parameter')
    p.add_argument('--N', type=int, default=8, help='number of frequency modes')
    p.add_argument('--steps', type=int, default=51, help='number of t values to sample')
    p.add_argument('--plot', action='store_true', help='generate plots')
    args = p.parse_args()

    # Run the trajectory tracking
    result = track_eigenvalue_trajectory(args.c, args.N, num_steps=args.steps)

    # Analyze for avoided crossings
    crossings = analyze_avoided_crossings(result, tol=1e-2)

    # Extract scaling laws
    extract_scaling_laws_from_trajectory(result)

    # Save data
    output_json = f'homotopy_trajectory_c{args.c}_N{args.N}.json'
    with open(output_json, 'w') as f:
        data_to_save = {
            'c': args.c,
            'N': args.N,
            't_values': result['t_values'].tolist(),
            'trajectories': {k: (v[0].tolist(), v[1].tolist()) for k, v in result['trajectories'].items()},
            'crossings': crossings,
        }
        json.dump(data_to_save, f, indent=2)
    print(f"\nSaved trajectory data to {output_json}")

    # Plot if requested
    if args.plot:
        plot_trajectories(result)

    print("\n" + "="*100)
    print("HOMOTOPY TRAJECTORY ANALYSIS COMPLETE")
    print("="*100)


if __name__ == '__main__':
    main()

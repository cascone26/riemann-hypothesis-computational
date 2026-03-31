"""
Fine-Grained Refinement Search Around Winning Operator

Winner from parameterized search:
  H = U(x)p + V(x)/p
  U(x) = 4.63 * x^0.876
  V(x) = 19.16 * x^(-0.167)
  Score: 0.004 (first 50 zeros)

This script does a systematic fine-grained search around those parameters,
using the first 100 zeros and optimizing (a, b) at each (alpha, beta) grid point.

Key question: Is alpha=0.876 genuinely special, or does alpha=1.0
(Berry-Keating) work just as well with the right coefficients?
"""

import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


# =============================================================================
# LOAD ZEROS
# =============================================================================

def load_zeros(count=100):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


# =============================================================================
# WKB COUNTING FUNCTION
# =============================================================================

def N_wkb(E, a, alpha, b, beta, x_range=(0.01, 300)):
    """
    WKB counting function for H = U(x)p + V(x)/p
    where U(x) = a*x^alpha, V(x) = b*x^beta.

    N_WKB(E) = (1/2pi) * integral of sqrt(E^2 - 4*U(x)*V(x)) / U(x) dx
    over the classically allowed region where E^2 >= 4*U(x)*V(x).
    """
    def integrand(x):
        U = a * x**alpha
        V = b * x**beta
        disc = E**2 - 4 * U * V
        if disc <= 0 or U <= 0:
            return 0.0
        return np.sqrt(disc) / U

    try:
        result, _ = quad(integrand, x_range[0], x_range[1],
                         limit=200, epsabs=1e-8, epsrel=1e-8)
        return result / (2 * np.pi)
    except:
        return 0.0


def N_wkb_vectorized(zeros, a, alpha, b, beta):
    """Compute N_WKB at each zero location."""
    return np.array([N_wkb(E, a, alpha, b, beta) for E in zeros])


# =============================================================================
# SCORING
# =============================================================================

def score(a, alpha, b, beta, zeros):
    """
    Score = (1/n) * sum((N_WKB(t_n) - n)^2)
    where t_n is the n-th zero and we expect N_WKB(t_n) ~ n.
    """
    n_zeros = len(zeros)
    total = 0.0
    for i, t_n in enumerate(zeros):
        n = i + 1  # expected count at the n-th zero
        n_wkb = N_wkb(t_n, a, alpha, b, beta)
        total += (n_wkb - n)**2
    return total / n_zeros


def score_ab(params, alpha, beta, zeros):
    """Score function for optimizing (a, b) at fixed (alpha, beta)."""
    a, b = params
    if a <= 0 or b <= 0:
        return 1e10
    return score(a, alpha, b, beta, zeros)


# =============================================================================
# GRID SEARCH WITH OPTIMIZATION
# =============================================================================

def refine_search(zeros, n_alpha=50, n_beta=50):
    """
    Fine-grained search over (alpha, beta) grid.
    At each grid point, optimize (a, b) using Nelder-Mead.
    """
    alphas = np.linspace(0.7, 1.0, n_alpha)
    betas = np.linspace(-0.5, 0.5, n_beta)

    score_grid = np.full((n_alpha, n_beta), np.nan)
    best_ab_grid = np.full((n_alpha, n_beta, 2), np.nan)

    global_best_score = float('inf')
    global_best_params = None

    total = n_alpha * n_beta
    start_time = time.time()

    for i, alpha in enumerate(alphas):
        for j, beta in enumerate(betas):
            # Starting guess near the known winner
            a0, b0 = 4.63, 19.16

            try:
                opt = minimize(score_ab, [a0, b0],
                               args=(alpha, beta, zeros),
                               method='Nelder-Mead',
                               options={'maxiter': 300, 'xatol': 1e-4, 'fatol': 1e-6})
                s = opt.fun
                a_opt, b_opt = opt.x
            except:
                s = 1e10
                a_opt, b_opt = a0, b0

            # Also try a wider starting point
            try:
                opt2 = minimize(score_ab, [1.0, 1.0],
                                args=(alpha, beta, zeros),
                                method='Nelder-Mead',
                                options={'maxiter': 300, 'xatol': 1e-4, 'fatol': 1e-6})
                if opt2.fun < s:
                    s = opt2.fun
                    a_opt, b_opt = opt2.x
            except:
                pass

            score_grid[i, j] = s
            best_ab_grid[i, j, 0] = a_opt
            best_ab_grid[i, j, 1] = b_opt

            if s < global_best_score:
                global_best_score = s
                global_best_params = (alpha, beta, a_opt, b_opt, s)

            done = i * n_beta + j + 1
            if done % 100 == 0 or done == total:
                elapsed = time.time() - start_time
                rate = done / elapsed if elapsed > 0 else 0
                eta = (total - done) / rate if rate > 0 else 0
                print(f"  [{done}/{total}] best={global_best_score:.6f} "
                      f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")

    return alphas, betas, score_grid, best_ab_grid, global_best_params


# =============================================================================
# BERRY-KEATING COMPARISON
# =============================================================================

def check_berry_keating(zeros):
    """
    Check alpha=1.0 specifically (Berry-Keating limit).
    Optimize over (beta, a, b).
    """
    print("\n" + "=" * 60)
    print("BERRY-KEATING CHECK: alpha = 1.0")
    print("=" * 60)

    betas_fine = np.linspace(-0.5, 0.5, 100)
    results = []

    for beta in betas_fine:
        try:
            opt = minimize(score_ab, [4.63, 19.16],
                           args=(1.0, beta, zeros),
                           method='Nelder-Mead',
                           options={'maxiter': 300})
            results.append((opt.fun, beta, opt.x[0], opt.x[1]))
        except:
            results.append((1e10, beta, 0, 0))

    results.sort(key=lambda x: x[0])
    best = results[0]
    print(f"  Best at alpha=1.0: beta={best[1]:.4f}, a={best[2]:.4f}, b={best[3]:.4f}")
    print(f"  Score: {best[0]:.6f}")
    return best


# =============================================================================
# PLOTTING
# =============================================================================

def plot_heatmap(alphas, betas, score_grid, global_best, bk_result, zeros):
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # --- Panel 1: Heatmap ---
    ax = axes[0, 0]
    # Use log scale for better visualization
    log_scores = np.log10(np.clip(score_grid, 1e-6, None))
    im = ax.pcolormesh(betas, alphas, log_scores, cmap='viridis_r', shading='auto')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('log10(score)')
    ax.set_xlabel('beta')
    ax.set_ylabel('alpha')
    ax.set_title('Score vs (alpha, beta) — optimized (a, b) at each point')

    # Mark the winner and Berry-Keating
    alpha_best, beta_best = global_best[0], global_best[1]
    ax.plot(beta_best, alpha_best, 'r*', markersize=15, label=f'Best: ({alpha_best:.3f}, {beta_best:.3f})')
    ax.axhline(y=1.0, color='white', linestyle='--', alpha=0.5, label='alpha=1 (Berry-Keating)')
    ax.axhline(y=0.876, color='orange', linestyle='--', alpha=0.5, label='alpha=0.876 (original)')
    ax.legend(fontsize=8, loc='lower right')

    # --- Panel 2: Score vs alpha (marginal, best over beta) ---
    ax = axes[0, 1]
    best_over_beta = np.nanmin(score_grid, axis=1)
    ax.semilogy(alphas, best_over_beta, 'b-', linewidth=2)
    ax.axvline(x=0.876, color='orange', linestyle='--', alpha=0.7, label='alpha=0.876')
    ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.7, label='alpha=1.0 (BK)')
    ax.axvline(x=alpha_best, color='green', linestyle='--', alpha=0.7, label=f'best={alpha_best:.3f}')
    ax.set_xlabel('alpha')
    ax.set_ylabel('Best score (log scale)')
    ax.set_title('Best score vs alpha (minimized over beta, a, b)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # --- Panel 3: Zero-by-zero comparison ---
    ax = axes[1, 0]
    alpha_b, beta_b, a_b, b_b = global_best[:4]
    n_show = 30
    ns = np.arange(1, n_show + 1)
    n_wkb_vals = N_wkb_vectorized(zeros[:n_show], a_b, alpha_b, b_b, beta_b)

    ax.plot(ns, ns, 'k--', linewidth=1, label='Perfect: N=n')
    ax.plot(ns, n_wkb_vals, 'bo-', markersize=4, linewidth=1, label='N_WKB(t_n)')
    ax.set_xlabel('Zero index n')
    ax.set_ylabel('N_WKB(t_n)')
    ax.set_title(f'Zero-by-zero: best params (a={a_b:.2f}, alpha={alpha_b:.3f}, b={b_b:.2f}, beta={beta_b:.3f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Panel 4: Residuals ---
    ax = axes[1, 1]
    residuals = n_wkb_vals - ns
    ax.bar(ns, residuals, color='steelblue', alpha=0.7)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.set_xlabel('Zero index n')
    ax.set_ylabel('N_WKB(t_n) - n')
    ax.set_title(f'Residuals (score={global_best[4]:.6f})')
    ax.grid(True, alpha=0.3)

    # Also show Berry-Keating residuals
    bk_s, bk_beta, bk_a, bk_b = bk_result
    n_wkb_bk = N_wkb_vectorized(zeros[:n_show], bk_a, 1.0, bk_b, bk_beta)
    residuals_bk = n_wkb_bk - ns
    ax.bar(ns + 0.3, residuals_bk, width=0.3, color='red', alpha=0.5, label=f'BK alpha=1 (score={bk_s:.4f})')
    ax.legend(fontsize=8)

    plt.suptitle('Refinement Search: H = a*x^alpha * p + b*x^beta / p',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'search_refinement.png'), dpi=150)
    plt.close()
    print(f"\nPlot saved to {os.path.join(PLOTS_DIR, 'search_refinement.png')}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("FINE-GRAINED REFINEMENT SEARCH")
    print("=" * 60)
    print("Searching around winning operator: U=4.63*x^0.876, V=19.16*x^(-0.167)")
    print("Using first 100 zeros for scoring\n")

    zeros = load_zeros(100)
    print(f"Loaded {len(zeros)} zeros (range: {zeros[0]:.4f} to {zeros[-1]:.4f})")

    # --- Main grid search ---
    print(f"\nGrid search: alpha in [0.7, 1.0] x beta in [-0.5, 0.5], 50x50 = 2500 points")
    print("At each point, optimizing (a, b) via Nelder-Mead...\n")

    alphas, betas, score_grid, best_ab_grid, global_best = refine_search(zeros)
    alpha_best, beta_best, a_best, b_best, s_best = global_best

    print(f"\n{'=' * 60}")
    print(f"GLOBAL BEST:")
    print(f"  alpha = {alpha_best:.6f}")
    print(f"  beta  = {beta_best:.6f}")
    print(f"  a     = {a_best:.6f}")
    print(f"  b     = {b_best:.6f}")
    print(f"  score = {s_best:.6f}")
    print(f"{'=' * 60}")

    # --- Berry-Keating comparison ---
    bk_result = check_berry_keating(zeros)
    bk_score = bk_result[0]

    print(f"\n{'=' * 60}")
    print(f"COMPARISON: alpha=0.876 vs alpha=1.0 (Berry-Keating)")
    print(f"{'=' * 60}")

    # Score at exact alpha=0.876
    opt_876 = minimize(score_ab, [4.63, 19.16],
                       args=(0.876, -0.167, zeros),
                       method='Nelder-Mead',
                       options={'maxiter': 500})
    s_876 = opt_876.fun
    a_876, b_876 = opt_876.x

    print(f"  alpha=0.876, beta=-0.167 (original): score={s_876:.6f} (a={a_876:.4f}, b={b_876:.4f})")
    print(f"  alpha=1.0   (Berry-Keating best):     score={bk_score:.6f}")
    print(f"  Global best:                           score={s_best:.6f}")
    ratio = bk_score / s_best if s_best > 0 else float('inf')
    print(f"\n  Berry-Keating is {ratio:.1f}x worse than global best")
    if ratio > 2:
        print(f"  --> alpha != 1 is SIGNIFICANTLY better. Berry-Keating is NOT optimal.")
    elif ratio > 1.2:
        print(f"  --> alpha != 1 is moderately better.")
    else:
        print(f"  --> Berry-Keating is competitive. alpha=1 may be sufficient.")

    # --- Zero-by-zero comparison for first 30 zeros ---
    print(f"\n{'=' * 60}")
    print(f"ZERO-BY-ZERO COMPARISON (first 30 zeros)")
    print(f"{'=' * 60}")
    print(f"  {'n':>3} {'t_n':>12} {'N_WKB(best)':>12} {'N_WKB(BK)':>12} {'err_best':>10} {'err_BK':>10}")
    print(f"  {'-' * 63}")

    bk_s, bk_beta, bk_a, bk_b = bk_result
    for i in range(30):
        t_n = zeros[i]
        n = i + 1
        n_best = N_wkb(t_n, a_best, alpha_best, b_best, beta_best)
        n_bk = N_wkb(t_n, bk_a, 1.0, bk_b, bk_beta)
        err_best = n_best - n
        err_bk = n_bk - n
        print(f"  {n:>3} {t_n:>12.4f} {n_best:>12.4f} {n_bk:>12.4f} {err_best:>+10.4f} {err_bk:>+10.4f}")

    # --- Check: how special is the valley in alpha? ---
    print(f"\n{'=' * 60}")
    print(f"ALPHA SENSITIVITY (best score at each alpha, optimized over everything else)")
    print(f"{'=' * 60}")

    alpha_check = [0.70, 0.75, 0.80, 0.85, 0.876, 0.90, 0.95, 1.00]
    for ac in alpha_check:
        # Find best beta for this alpha from the grid
        idx = np.argmin(np.abs(alphas - ac))
        best_j = np.nanargmin(score_grid[idx, :])
        s_here = score_grid[idx, best_j]
        a_here = best_ab_grid[idx, best_j, 0]
        b_here = best_ab_grid[idx, best_j, 1]
        beta_here = betas[best_j]
        print(f"  alpha={ac:.3f}: score={s_here:.6f} (beta={beta_here:.3f}, a={a_here:.2f}, b={b_here:.2f})")

    # --- Plot ---
    plot_heatmap(alphas, betas, score_grid, global_best, bk_result, zeros)

    # --- Save results ---
    results = {
        "global_best": {
            "alpha": float(alpha_best),
            "beta": float(beta_best),
            "a": float(a_best),
            "b": float(b_best),
            "score": float(s_best),
        },
        "berry_keating": {
            "alpha": 1.0,
            "beta": float(bk_result[1]),
            "a": float(bk_result[2]),
            "b": float(bk_result[3]),
            "score": float(bk_result[0]),
        },
        "original_winner": {
            "alpha": 0.876,
            "beta": -0.167,
            "a": float(a_876),
            "b": float(b_876),
            "score": float(s_876),
        },
        "n_zeros_used": 100,
        "grid_size": "50x50",
    }

    import json
    with open(os.path.join(RESULTS_DIR, "refine_search.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {os.path.join(RESULTS_DIR, 'refine_search.json')}")


if __name__ == "__main__":
    main()

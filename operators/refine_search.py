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

Performance: Uses numpy trapezoid integration instead of scipy.quad for speed.
"""

import numpy as np
from scipy.optimize import minimize
import os
import sys
import time
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Precompute integration grid (log-spaced for better resolution near origin)
X_GRID = np.logspace(np.log10(0.01), np.log10(500), 2000)
X_DIFF = np.diff(X_GRID)  # for trapezoid rule


def flush_print(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()


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
# FAST WKB COUNTING FUNCTION (vectorized numpy)
# =============================================================================

def N_wkb_fast(E, a, alpha, b, beta):
    """
    WKB counting function for H = U(x)p + V(x)/p
    where U(x) = a*x^alpha, V(x) = b*x^beta.

    N_WKB(E) = (1/2pi) * integral of sqrt(E^2 - 4*U(x)*V(x)) / U(x) dx
    over the classically allowed region where E^2 >= 4*U(x)*V(x).

    Uses precomputed X_GRID with trapezoid rule for speed.
    """
    U = a * X_GRID**alpha
    V = b * X_GRID**beta
    disc = E**2 - 4 * U * V
    # Classically allowed region: disc > 0
    mask = disc > 0
    integrand = np.zeros_like(X_GRID)
    integrand[mask] = np.sqrt(disc[mask]) / U[mask]
    # Trapezoid rule
    result = np.trapezoid(integrand, X_GRID)
    return result / (2 * np.pi)


def N_wkb_batch(zeros, a, alpha, b, beta):
    """Compute N_WKB at all zero locations. Returns array."""
    # Precompute U, V on the grid once
    U = a * X_GRID**alpha
    V = b * X_GRID**beta
    product_4UV = 4 * U * V

    results = np.empty(len(zeros))
    for i, E in enumerate(zeros):
        disc = E**2 - product_4UV
        mask = disc > 0
        integrand = np.zeros_like(X_GRID)
        integrand[mask] = np.sqrt(disc[mask]) / U[mask]
        results[i] = np.trapezoid(integrand, X_GRID) / (2 * np.pi)
    return results


# =============================================================================
# SCORING
# =============================================================================

def score_fast(a, alpha, b, beta, zeros):
    """
    Score = (1/n) * sum((N_WKB(t_n) - n)^2)
    where t_n is the n-th zero and we expect N_WKB(t_n) ~ n.
    """
    if a <= 0 or b <= 0:
        return 1e10
    try:
        n_wkb = N_wkb_batch(zeros, a, alpha, b, beta)
        ns = np.arange(1, len(zeros) + 1, dtype=float)
        return np.mean((n_wkb - ns)**2)
    except:
        return 1e10


def score_ab(params, alpha, beta, zeros):
    """Score function for optimizing (a, b) at fixed (alpha, beta)."""
    a, b = params
    return score_fast(a, alpha, b, beta, zeros)


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
                               options={'maxiter': 200, 'xatol': 1e-3, 'fatol': 1e-5})
                s = opt.fun
                a_opt, b_opt = opt.x
            except:
                s = 1e10
                a_opt, b_opt = a0, b0

            # Also try a second starting point
            try:
                opt2 = minimize(score_ab, [1.0, 5.0],
                                args=(alpha, beta, zeros),
                                method='Nelder-Mead',
                                options={'maxiter': 200, 'xatol': 1e-3, 'fatol': 1e-5})
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
            if done % 50 == 0 or done == total:
                elapsed = time.time() - start_time
                rate = done / elapsed if elapsed > 0 else 0
                eta = (total - done) / rate if rate > 0 else 0
                flush_print(f"  [{done}/{total}] best={global_best_score:.6f} "
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
    flush_print("\n" + "=" * 60)
    flush_print("BERRY-KEATING CHECK: alpha = 1.0")
    flush_print("=" * 60)

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
    flush_print(f"  Best at alpha=1.0: beta={best[1]:.4f}, a={best[2]:.4f}, b={best[3]:.4f}")
    flush_print(f"  Score: {best[0]:.6f}")
    return best


# =============================================================================
# PLOTTING
# =============================================================================

def plot_heatmap(alphas, betas, score_grid, global_best, bk_result, zeros):
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # --- Panel 1: Heatmap ---
    ax = axes[0, 0]
    log_scores = np.log10(np.clip(score_grid, 1e-6, None))
    im = ax.pcolormesh(betas, alphas, log_scores, cmap='viridis_r', shading='auto')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('log10(score)')
    ax.set_xlabel('beta')
    ax.set_ylabel('alpha')
    ax.set_title('Score vs (alpha, beta) -- optimized (a, b) at each point')

    alpha_best, beta_best = global_best[0], global_best[1]
    ax.plot(beta_best, alpha_best, 'r*', markersize=15,
            label=f'Best: ({alpha_best:.3f}, {beta_best:.3f})')
    ax.axhline(y=1.0, color='white', linestyle='--', alpha=0.5, label='alpha=1 (Berry-Keating)')
    ax.axhline(y=0.876, color='orange', linestyle='--', alpha=0.5, label='alpha=0.876 (original)')
    ax.legend(fontsize=8, loc='lower right')

    # --- Panel 2: Score vs alpha (marginal, best over beta) ---
    ax = axes[0, 1]
    best_over_beta = np.nanmin(score_grid, axis=1)
    ax.semilogy(alphas, best_over_beta, 'b-', linewidth=2)
    ax.axvline(x=0.876, color='orange', linestyle='--', alpha=0.7, label='alpha=0.876')
    ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.7, label='alpha=1.0 (BK)')
    ax.axvline(x=alpha_best, color='green', linestyle='--', alpha=0.7,
               label=f'best={alpha_best:.3f}')
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
    n_wkb_vals = N_wkb_batch(zeros[:n_show], a_b, alpha_b, b_b, beta_b)

    ax.plot(ns, ns, 'k--', linewidth=1, label='Perfect: N=n')
    ax.plot(ns, n_wkb_vals, 'bo-', markersize=4, linewidth=1, label='N_WKB(t_n)')
    ax.set_xlabel('Zero index n')
    ax.set_ylabel('N_WKB(t_n)')
    ax.set_title(f'Zero-by-zero: a={a_b:.2f}, alpha={alpha_b:.3f}, '
                 f'b={b_b:.2f}, beta={beta_b:.3f}')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Panel 4: Residuals ---
    ax = axes[1, 1]
    residuals = n_wkb_vals - ns
    ax.bar(ns - 0.15, residuals, width=0.3, color='steelblue', alpha=0.7,
           label=f'Best (score={global_best[4]:.4f})')
    ax.axhline(y=0, color='k', linewidth=0.5)

    # Berry-Keating residuals
    bk_s, bk_beta, bk_a, bk_b = bk_result
    n_wkb_bk = N_wkb_batch(zeros[:n_show], bk_a, 1.0, bk_b, bk_beta)
    residuals_bk = n_wkb_bk - ns
    ax.bar(ns + 0.15, residuals_bk, width=0.3, color='red', alpha=0.5,
           label=f'BK alpha=1 (score={bk_s:.4f})')

    ax.set_xlabel('Zero index n')
    ax.set_ylabel('N_WKB(t_n) - n')
    ax.set_title('Residuals: Best vs Berry-Keating')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.suptitle('Refinement Search: H = a*x^alpha * p + b*x^beta / p',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'search_refinement.png'), dpi=150)
    plt.close()
    flush_print(f"\nPlot saved to {os.path.join(PLOTS_DIR, 'search_refinement.png')}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    flush_print("=" * 60)
    flush_print("FINE-GRAINED REFINEMENT SEARCH")
    flush_print("=" * 60)
    flush_print("Searching around winning operator: U=4.63*x^0.876, V=19.16*x^(-0.167)")
    flush_print("Using first 100 zeros for scoring\n")

    zeros = load_zeros(100)
    flush_print(f"Loaded {len(zeros)} zeros (range: {zeros[0]:.4f} to {zeros[-1]:.4f})")

    # Quick sanity check: score the original winner
    s_orig = score_fast(4.63, 0.876, 19.16, -0.167, zeros)
    flush_print(f"Original winner score (100 zeros): {s_orig:.6f}")

    # --- Main grid search ---
    flush_print(f"\nGrid search: alpha in [0.7, 1.0] x beta in [-0.5, 0.5], 50x50 = 2500 points")
    flush_print("At each point, optimizing (a, b) via Nelder-Mead...\n")

    alphas, betas, score_grid, best_ab_grid, global_best = refine_search(zeros)
    alpha_best, beta_best, a_best, b_best, s_best = global_best

    flush_print(f"\n{'=' * 60}")
    flush_print(f"GLOBAL BEST:")
    flush_print(f"  alpha = {alpha_best:.6f}")
    flush_print(f"  beta  = {beta_best:.6f}")
    flush_print(f"  a     = {a_best:.6f}")
    flush_print(f"  b     = {b_best:.6f}")
    flush_print(f"  score = {s_best:.6f}")
    flush_print(f"{'=' * 60}")

    # --- Berry-Keating comparison ---
    bk_result = check_berry_keating(zeros)
    bk_score = bk_result[0]

    flush_print(f"\n{'=' * 60}")
    flush_print(f"COMPARISON: alpha=0.876 vs alpha=1.0 (Berry-Keating)")
    flush_print(f"{'=' * 60}")

    # Score at exact alpha=0.876 with optimized a,b
    opt_876 = minimize(score_ab, [4.63, 19.16],
                       args=(0.876, -0.167, zeros),
                       method='Nelder-Mead',
                       options={'maxiter': 500})
    s_876 = opt_876.fun
    a_876, b_876 = opt_876.x

    flush_print(f"  alpha=0.876, beta=-0.167 (original): score={s_876:.6f} "
                f"(a={a_876:.4f}, b={b_876:.4f})")
    flush_print(f"  alpha=1.0   (Berry-Keating best):     score={bk_score:.6f}")
    flush_print(f"  Global best:                           score={s_best:.6f}")
    ratio = bk_score / s_best if s_best > 0 else float('inf')
    flush_print(f"\n  Berry-Keating is {ratio:.1f}x worse than global best")
    if ratio > 2:
        flush_print("  --> alpha != 1 is SIGNIFICANTLY better. Berry-Keating is NOT optimal.")
    elif ratio > 1.2:
        flush_print("  --> alpha != 1 is moderately better.")
    else:
        flush_print("  --> Berry-Keating is competitive. alpha=1 may be sufficient.")

    # --- Zero-by-zero comparison for first 30 zeros ---
    flush_print(f"\n{'=' * 60}")
    flush_print(f"ZERO-BY-ZERO COMPARISON (first 30 zeros)")
    flush_print(f"{'=' * 60}")
    flush_print(f"  {'n':>3} {'t_n':>12} {'N_WKB(best)':>12} {'N_WKB(BK)':>12} "
                f"{'err_best':>10} {'err_BK':>10}")
    flush_print(f"  {'-' * 63}")

    bk_s, bk_beta, bk_a, bk_b = bk_result
    n_wkb_best_30 = N_wkb_batch(zeros[:30], a_best, alpha_best, b_best, beta_best)
    n_wkb_bk_30 = N_wkb_batch(zeros[:30], bk_a, 1.0, bk_b, bk_beta)

    for i in range(30):
        t_n = zeros[i]
        n = i + 1
        n_best = n_wkb_best_30[i]
        n_bk = n_wkb_bk_30[i]
        err_best = n_best - n
        err_bk = n_bk - n
        flush_print(f"  {n:>3} {t_n:>12.4f} {n_best:>12.4f} {n_bk:>12.4f} "
                     f"{err_best:>+10.4f} {err_bk:>+10.4f}")

    # --- Check: how special is the valley in alpha? ---
    flush_print(f"\n{'=' * 60}")
    flush_print(f"ALPHA SENSITIVITY (best score at each alpha, optimized over everything else)")
    flush_print(f"{'=' * 60}")

    alpha_check = [0.70, 0.75, 0.80, 0.85, 0.876, 0.90, 0.95, 1.00]
    for ac in alpha_check:
        idx = np.argmin(np.abs(alphas - ac))
        best_j = np.nanargmin(score_grid[idx, :])
        s_here = score_grid[idx, best_j]
        a_here = best_ab_grid[idx, best_j, 0]
        b_here = best_ab_grid[idx, best_j, 1]
        beta_here = betas[best_j]
        flush_print(f"  alpha={ac:.3f}: score={s_here:.6f} "
                     f"(beta={beta_here:.3f}, a={a_here:.2f}, b={b_here:.2f})")

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

    with open(os.path.join(RESULTS_DIR, "refine_search.json"), "w") as f:
        json.dump(results, f, indent=2)
    flush_print(f"\nResults saved to {os.path.join(RESULTS_DIR, 'refine_search.json')}")


if __name__ == "__main__":
    main()

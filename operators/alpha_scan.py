"""
Scan alpha to understand what makes alpha ≈ 0.876 special.

For each alpha, optimize (a, b, beta) and compute the best score.
Plot score vs alpha to see if there's a sharp minimum.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")


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


def N_wkb(E, alpha, beta, a, b):
    def integrand(x):
        U = a * x**alpha
        V = b * x**beta
        disc = E**2 - 4*U*V
        if disc <= 0 or U <= 0:
            return 0.0
        return math.sqrt(disc) / U
    try:
        result, _ = quad(integrand, 0.01, 200, limit=100, epsabs=1e-6)
        return result / (2 * math.pi)
    except:
        return 0.0


def score(alpha, beta, a, b, zeros):
    total = 0.0
    for i, E in enumerate(zeros[:50]):
        n_wkb = N_wkb(E, alpha, beta, a, b)
        n_exact = i + 1
        total += (n_wkb - n_exact)**2
    return total / 50


def optimize_for_alpha(alpha, zeros):
    """For fixed alpha, find best (beta, a, b)."""
    def objective(params):
        beta, a, b = params
        if a <= 0 or b <= 0:
            return 1e10
        return score(alpha, beta, a, b, zeros)

    best = float('inf')
    best_params = None

    # Try several starting points
    for _ in range(20):
        beta0 = np.random.uniform(-0.5, 0.5)
        a0 = np.random.uniform(0.5, 10)
        b0 = np.random.uniform(0.5, 30)
        try:
            result = minimize(objective, [beta0, a0, b0],
                            method='Nelder-Mead',
                            options={'maxiter': 300})
            if result.fun < best:
                best = result.fun
                best_params = result.x
        except:
            pass

    return best, best_params


def main():
    print("ALPHA DEPENDENCE SCAN")
    print("=" * 50)

    zeros = load_zeros(100)
    np.random.seed(42)

    alphas = np.linspace(0.5, 1.5, 40)
    scores_list = []
    params_list = []

    for alpha in alphas:
        s, params = optimize_for_alpha(alpha, zeros)
        scores_list.append(s)
        params_list.append(params)
        beta, a, b = params if params is not None else (0, 0, 0)
        print(f"  alpha={alpha:.3f}: score={s:.6f}, beta={beta:.3f}, a={a:.3f}, b={b:.3f}")

    scores_arr = np.array(scores_list)
    best_idx = np.argmin(scores_arr)
    best_alpha = alphas[best_idx]

    print(f"\nBest alpha: {best_alpha:.4f} (score: {scores_arr[best_idx]:.6f})")
    print(f"Berry-Keating (alpha=1.0): score = {scores_arr[np.argmin(np.abs(alphas - 1.0))]:.6f}")

    # Is alpha=1 a local minimum?
    idx_1 = np.argmin(np.abs(alphas - 1.0))
    if idx_1 > 0 and idx_1 < len(alphas) - 1:
        local_min = scores_arr[idx_1] < scores_arr[idx_1-1] and scores_arr[idx_1] < scores_arr[idx_1+1]
        print(f"alpha=1.0 is a local minimum: {local_min}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.plot(alphas, scores_arr, 'b.-', markersize=4)
    ax.axvline(x=best_alpha, color='red', linestyle='--', alpha=0.5,
               label=f'Best: alpha={best_alpha:.3f}')
    ax.axvline(x=1.0, color='green', linestyle=':', alpha=0.5,
               label='Berry-Keating: alpha=1.0')
    ax.set_xlabel('alpha')
    ax.set_ylabel('Score (lower = better)')
    ax.set_title('Score vs Alpha (other params optimized)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')

    ax = axes[1]
    betas = [p[0] if p is not None else 0 for p in params_list]
    ax.plot(alphas, betas, 'r.-', markersize=4)
    ax.axvline(x=best_alpha, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('alpha')
    ax.set_ylabel('Optimal beta')
    ax.set_title('Optimal Beta vs Alpha')
    ax.grid(True, alpha=0.3)

    plt.suptitle('Alpha Dependence — What Makes alpha≈0.876 Special?', fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'alpha_scan.png'), dpi=150)
    plt.close()
    print(f"\nPlot saved.")

    # Save
    import json
    results = {
        "best_alpha": float(best_alpha),
        "best_score": float(scores_arr[best_idx]),
        "alpha_1_score": float(scores_arr[idx_1]),
        "improvement_over_BK": float(scores_arr[idx_1] / scores_arr[best_idx]),
        "scan": [{"alpha": float(a), "score": float(s)} for a, s in zip(alphas, scores_arr)],
    }
    with open(os.path.join(RESULTS_DIR, "alpha_scan.json"), "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()

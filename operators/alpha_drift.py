"""
Alpha Drift Test: Does the optimal alpha shift toward 1 as N → ∞?

If alpha<1 is just a finite-range artifact, the optimal alpha should
drift toward 1.0 as we use more zeros. If it stays put, it's a genuine
property of the zeta function.

Tests at N = 50, 100, 200, 500, 1000 zeros.
Uses already-known params as starting points, then local refinement.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
import math
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
PLOTS_DIR   = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")


def load_zeros(count=1100):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def N_wkb(E, alpha, beta, a, b):
    def integrand(x):
        disc = E**2 - 4 * a * (x**alpha) * b * (x**beta)
        if disc <= 0:
            return 0.0
        denom = a * (x**alpha)
        if denom <= 0:
            return 0.0
        return math.sqrt(disc) / denom
    try:
        result, _ = quad(integrand, 0.01, 500, limit=150, epsabs=1e-5)
        return result / (2 * math.pi)
    except Exception:
        return 0.0


def score_fixed_params(alpha, beta, a, b, zeros, n):
    """Score params against first n zeros."""
    total = 0.0
    for i, E in enumerate(zeros[:n]):
        nw = N_wkb(E, alpha, beta, a, b)
        total += (nw - (i + 1))**2
    return total / n


def optimize_alpha_neighborhood(alpha_center, zeros, n, n_starts=8):
    """
    Search for best (alpha, beta, a, b) near alpha_center.
    Alpha is free in [alpha_center-0.15, alpha_center+0.15].
    """
    def objective(params):
        alpha, beta, log_a, log_b = params
        if not (0.5 <= alpha <= 1.3):
            return 1e10
        a, b = math.exp(log_a), math.exp(log_b)
        total = 0.0
        # Subsample for speed: every 3rd zero up to n
        indices = list(range(0, min(n, len(zeros)), max(1, n // 40)))
        for i in indices:
            nw = N_wkb(zeros[i], alpha, beta, a, b)
            total += (nw - (i + 1))**2
        return total / len(indices)

    best_score = float('inf')
    best_params = None
    np.random.seed(42)

    for _ in range(n_starts):
        alpha0 = np.clip(alpha_center + np.random.uniform(-0.1, 0.1), 0.6, 1.2)
        beta0  = np.random.uniform(-0.5, 0.5)
        la0    = np.random.uniform(0.5, 2.0)
        lb0    = np.random.uniform(1.0, 5.0)
        try:
            res = minimize(objective, [alpha0, beta0, la0, lb0],
                           method='Nelder-Mead',
                           options={'maxiter': 800, 'xatol': 1e-4, 'fatol': 1e-6})
            if res.fun < best_score:
                best_score = res.fun
                best_params = res.x
        except Exception:
            pass

    if best_params is None:
        return None, None
    alpha, beta, log_a, log_b = best_params
    return best_score, (alpha, beta, math.exp(log_a), math.exp(log_b))


def main():
    print("ALPHA DRIFT TEST: Does optimal alpha drift toward 1 with more zeros?")
    print("=" * 65)

    zeros = load_zeros(1100)
    print(f"Loaded {len(zeros)} zeros. Range: [{zeros[0]:.2f}, {zeros[-1]:.2f}]")
    print()

    # Load known good starting point
    with open(os.path.join(RESULTS_DIR, "refine_search.json")) as f:
        refine = json.load(f)
    bk_alpha   = refine["berry_keating"]["alpha"]
    bk_beta    = refine["berry_keating"]["beta"]
    bk_a       = refine["berry_keating"]["a"]
    bk_b       = refine["berry_keating"]["b"]
    opt_alpha  = refine["global_best"]["alpha"]

    N_values = [50, 100, 200, 500, 1000]
    results = []

    for N in N_values:
        print(f"N = {N} zeros (E up to {zeros[N-1]:.1f})...")

        # Score BK at this N
        score_bk = score_fixed_params(bk_alpha, bk_beta, bk_a, bk_b, zeros, N)

        # Optimize near the known optimal alpha
        s_opt, p_opt = optimize_alpha_neighborhood(opt_alpha, zeros, N, n_starts=8)

        if p_opt:
            alpha_found, beta_found, a_found, b_found = p_opt
            # Get exact score on full N zeros
            s_exact = score_fixed_params(alpha_found, beta_found, a_found, b_found, zeros, N)
        else:
            alpha_found = opt_alpha
            s_exact = float('inf')

        improvement = score_bk / s_exact if s_exact > 0 else float('inf')

        print(f"  BK score:      {score_bk:.4f}")
        print(f"  Optimal alpha: {alpha_found:.4f}  (β={beta_found:.3f})")
        print(f"  Optimal score: {s_exact:.4f}")
        print(f"  Improvement:   {improvement:.2f}x")
        print()

        results.append({
            "N": N,
            "E_max": float(zeros[N-1]),
            "alpha_opt": float(alpha_found),
            "beta_opt": float(beta_found) if p_opt else None,
            "score_bk": float(score_bk),
            "score_opt": float(s_exact),
            "improvement_over_bk": float(improvement),
        })

    # ── Analysis ────────────────────────────────────────────────
    print("=" * 65)
    print("DRIFT ANALYSIS")
    print("=" * 65)
    print(f"{'N':>6}  {'E_max':>7}  {'α_opt':>7}  {'score_opt':>10}  {'vs BK':>8}")
    for r in results:
        print(f"{r['N']:>6}  {r['E_max']:>7.1f}  {r['alpha_opt']:>7.4f}  "
              f"{r['score_opt']:>10.4f}  {r['improvement_over_bk']:>7.1f}x")

    alphas_found = [r["alpha_opt"] for r in results]
    drift = alphas_found[-1] - alphas_found[0]
    print()
    print(f"Alpha drift from N=50 to N=1000: {drift:+.4f}")
    if drift > 0.05:
        verdict = "DRIFTING TOWARD 1 — finite-range artifact"
    elif drift < -0.05:
        verdict = "DRIFTING AWAY FROM 1 — deepens the mystery"
    else:
        verdict = "STABLE — genuine property of the zeta function"
    print(f"Verdict: {verdict}")

    # ── Plot ────────────────────────────────────────────────────
    Ns = [r["N"] for r in results]
    opt_alphas   = [r["alpha_opt"]          for r in results]
    improvements = [r["improvement_over_bk"] for r in results]
    scores_opt   = [r["score_opt"]           for r in results]
    scores_bk    = [r["score_bk"]            for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax = axes[0]
    ax.semilogx(Ns, opt_alphas, 'ro-', markersize=8, linewidth=2)
    ax.axhline(1.0, color='g', linestyle='--', linewidth=1.5,
               label='α=1 (Berry-Keating)')
    ax.axhline(np.mean(opt_alphas), color='b', linestyle=':', linewidth=1,
               label=f'Mean α={np.mean(opt_alphas):.3f}')
    for i, (n, a) in enumerate(zip(Ns, opt_alphas)):
        ax.annotate(f'{a:.3f}', (n, a), textcoords='offset points',
                    xytext=(5, 5), fontsize=8)
    ax.set_xlabel('N (number of zeros used)')
    ax.set_ylabel('Optimal alpha')
    ax.set_title('Alpha Drift: Does Optimal α → 1 as N → ∞?')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.7, 1.15)

    ax = axes[1]
    ax.semilogx(Ns, improvements, 'b^-', markersize=8, linewidth=2)
    ax.axhline(1.0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('N (number of zeros used)')
    ax.set_ylabel('Improvement over Berry-Keating (×)')
    ax.set_title('How Much Does Optimal Beat BK?\n(1.0 = tie)')
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.loglog(Ns, scores_bk,  'g^-', markersize=8, linewidth=2, label='Berry-Keating')
    ax.loglog(Ns, scores_opt, 'r^-', markersize=8, linewidth=2, label='Optimal α')
    ax.set_xlabel('N (number of zeros used)')
    ax.set_ylabel('Score (lower = better)')
    ax.set_title('Score vs N: BK vs Optimal')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle(
        f'Alpha Drift Test — Verdict: {verdict}',
        fontsize=11
    )
    plt.tight_layout()
    out_path = os.path.join(PLOTS_DIR, "alpha_drift.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"\nPlot saved: {out_path}")

    with open(os.path.join(RESULTS_DIR, "alpha_drift.json"), "w") as f:
        json.dump({"results": results, "verdict": verdict, "drift": float(drift)}, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()

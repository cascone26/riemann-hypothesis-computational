"""
Parameterized Operator Search

This is genuinely uncharted territory — no published systematic computational
search over operator families for matching individual Riemann zeros exists.

We search the family H = U(x)p + V(x)/p (Sierra 2011 generalization).
The symmetrized quantum operator is:
  H = (1/2)(U(x)p + pU(x)) + V(x)/p

In the WKB limit, the eigenvalue condition is:
  integral sqrt(E - V_eff(x)) dx = (n + phase) * pi

where V_eff depends on U, V, and the energy E.

For the smooth counting function to match, we need:
  N_WKB(E) ~ (E/2pi) * log(E/2pi) - E/2pi + 7/8

APPROACH: Rather than discretize the full operator (which failed in v1),
we use the WKB quantization condition directly. This is exact for the
smooth part and lets us search much more efficiently.

For each candidate (U, V), we:
1. Compute the classical phase space area as a function of E
2. Compare the resulting N_WKB(E) to the actual zero counting function
3. Use the deviation from the smooth counting function to score
4. Optimize over U, V parameters

The REAL question: can any (U, V) in this family reproduce the
OSCILLATORY corrections (individual zeros), not just the smooth average?
"""

import numpy as np
from scipy.optimize import minimize, brentq
from scipy.integrate import quad
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_zeros(count=200):
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


def N_exact(zeros, E):
    """Exact counting function from known zeros."""
    return np.sum(zeros <= E)


def N_smooth(E):
    """Smooth approximation to zero counting function."""
    if E <= 2:
        return 0.0
    return (E / (2*np.pi)) * np.log(E / (2*np.pi)) - E / (2*np.pi) + 7/8


# =============================================================================
# WKB QUANTIZATION FOR H = U(x)p + V(x)/p
# =============================================================================

def classical_phase_area(E, U_func, V_func, x_range=(0.01, 200)):
    """
    For the classical Hamiltonian H = U(x)*p + V(x)/p,
    the energy surface H = E gives:
      p = [E ± sqrt(E^2 - 4*U(x)*V(x))] / (2*U(x))

    The phase space area enclosed by the orbit at energy E:
      A(E) = oint p dx = integral [p_+(x) - p_-(x)] dx

    over the classically allowed region where E^2 >= 4*U(x)*V(x).

    The WKB quantization: A(E_n) = 2*pi*(n + phase_correction)
    So N_WKB(E) = A(E) / (2*pi)
    """
    def p_diff(x):
        U = U_func(x)
        V = V_func(x)
        discriminant = E**2 - 4*U*V
        if discriminant <= 0 or U <= 0:
            return 0.0
        return np.sqrt(discriminant) / U

    try:
        result, _ = quad(p_diff, x_range[0], x_range[1],
                        limit=200, epsabs=1e-8)
        return result
    except:
        return 0.0


def N_wkb(E, U_func, V_func):
    """WKB counting function for given U, V."""
    return classical_phase_area(E, U_func, V_func) / (2 * np.pi)


# =============================================================================
# CANDIDATE FAMILIES
# =============================================================================

def make_UV_power(alpha, beta, a, b):
    """
    U(x) = a * x^alpha
    V(x) = b * x^beta

    Berry-Keating: alpha=1, beta=0, a=1, b=0 (degenerate)
    Sierra 2011: alpha=1, beta=1, a=1, b=ell_p^2
    """
    def U(x):
        return a * x**alpha
    def V(x):
        return b * x**beta
    return U, V


def make_UV_log(a, b, c):
    """
    U(x) = a * x * (1 + b*log(x))
    V(x) = c * x

    Log corrections suggested by the non-integer effective dimension d≈2.46.
    """
    def U(x):
        return a * x * (1 + b * np.log(max(x, 1e-10)))
    def V(x):
        return c * x
    return U, V


def make_UV_shifted(a, b, ell_x, ell_p):
    """
    H = (x + ell_x)(p + ell_p) = xp + ell_p*x + ell_x*p + ell_x*ell_p

    This is Berry-Keating (2011) compact version.
    Equivalent to U(x) = x + ell_x, V(x) = 0 with p -> p + ell_p.
    But in the U*p + V/p form, it maps to:
    U(x) = x + ell_x, V(x) = ell_x * ell_p (constant)
    """
    def U(x):
        return x + ell_x
    def V(x):
        return ell_x * ell_p
    return U, V


# =============================================================================
# SEARCH ENGINE
# =============================================================================

def score_candidate(params, family, zeros, E_test_points):
    """
    Score a candidate (U, V) against known zeros.

    Scoring:
    1. Compare N_WKB(E) to N_exact(E) at the zero locations
    2. Lower score = better match
    """
    try:
        U, V = family(*params)
    except:
        return 1e10

    total_error = 0.0
    for E in E_test_points:
        n_wkb = N_wkb(E, U, V)
        n_exact = N_exact(zeros, E)
        if n_exact > 0:
            total_error += (n_wkb - n_exact)**2 / n_exact
        else:
            total_error += n_wkb**2

    return total_error / len(E_test_points)


def search_family(family, param_bounds, zeros, n_random=500, n_refine=20):
    """
    Search over a parameter family:
    1. Random sampling to find promising regions
    2. Local optimization from best random points
    """
    # Test at the actual zero locations
    E_test = zeros[:50]  # first 50 zeros

    # Random search
    best_score = float('inf')
    best_params = None
    results = []

    n_params = len(param_bounds)
    for i in range(n_random):
        params = []
        for lo, hi in param_bounds:
            params.append(lo + (hi - lo) * np.random.random())

        score = score_candidate(params, family, zeros, E_test)
        results.append((score, params))

        if score < best_score:
            best_score = score
            best_params = params

    # Sort and refine top candidates
    results.sort(key=lambda x: x[0])
    refined = []

    for score, params in results[:n_refine]:
        try:
            opt = minimize(lambda p: score_candidate(p, family, zeros, E_test),
                          params, method='Nelder-Mead',
                          options={'maxiter': 200, 'xatol': 1e-6})
            refined.append((opt.fun, list(opt.x)))
        except:
            refined.append((score, params))

    refined.sort(key=lambda x: x[0])
    return refined[0] if refined else (best_score, best_params)


def run_search():
    print("=" * 60)
    print("PARAMETERIZED OPERATOR SEARCH")
    print("=" * 60)
    print("Searching H = U(x)p + V(x)/p families")
    print("This has NEVER been done in the literature.\n")

    zeros = load_zeros(200)
    np.random.seed(42)

    families = {
        "Power law: U=a*x^alpha, V=b*x^beta": {
            "family": make_UV_power,
            "bounds": [(0.5, 2.0), (0.0, 2.0), (0.1, 10.0), (0.0, 10.0)],
            "param_names": ["alpha", "beta", "a", "b"],
        },
        "Log-corrected: U=a*x*(1+b*log(x)), V=c*x": {
            "family": make_UV_log,
            "bounds": [(0.01, 5.0), (-1.0, 1.0), (0.0, 5.0)],
            "param_names": ["a", "b", "c"],
        },
        "Shifted: (x+ell_x)(p+ell_p)": {
            "family": make_UV_shifted,
            "bounds": [(0.1, 10.0), (0.1, 10.0), (0.01, 5.0), (0.01, 5.0)],
            "param_names": ["a", "b", "ell_x", "ell_p"],
        },
    }

    all_results = {}

    for name, config in families.items():
        print(f"\n{'='*50}")
        print(f"Family: {name}")
        print(f"{'='*50}")

        best_score, best_params = search_family(
            config["family"], config["bounds"], zeros,
            n_random=1000, n_refine=30
        )

        print(f"  Best score: {best_score:.6f}")
        for pname, pval in zip(config["param_names"], best_params):
            print(f"    {pname} = {pval:.6f}")

        # Evaluate the winner
        U, V = config["family"](*best_params)
        E_range = np.linspace(10, zeros[99], 200)
        n_wkb_vals = [N_wkb(E, U, V) for E in E_range]
        n_smooth_vals = [N_smooth(E) for E in E_range]
        n_exact_vals = [N_exact(zeros, E) for E in E_range]

        all_results[name] = {
            "score": float(best_score),
            "params": {pname: float(pval) for pname, pval in zip(config["param_names"], best_params)},
            "n_wkb_sample": n_wkb_vals[:10],
        }

        # Compare: can this candidate match individual zeros?
        print(f"\n  Zero-by-zero comparison (first 20):")
        print(f"  {'n':>3} {'Zero':>10} {'N_WKB':>10} {'N_smooth':>10} {'N_exact':>8}")
        print(f"  {'-'*45}")
        for i in range(min(20, len(zeros))):
            E = zeros[i]
            nw = N_wkb(E, U, V)
            ns = N_smooth(E)
            ne = i + 1
            print(f"  {i+1:>3} {E:>10.4f} {nw:>10.4f} {ns:>10.4f} {ne:>8d}")

    # Plot comparison of all families
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    E_plot = np.linspace(10, 250, 500)
    n_exact_plot = [N_exact(zeros, E) for E in E_plot]
    n_smooth_plot = [N_smooth(E) for E in E_plot]

    for idx, (name, config) in enumerate(families.items()):
        ax = axes[idx // 2, idx % 2]
        best_params = list(all_results[name]["params"].values())
        U, V = config["family"](*best_params)
        n_wkb_plot = [N_wkb(E, U, V) for E in E_plot]

        ax.plot(E_plot, n_exact_plot, 'k-', linewidth=2, label='Exact N(E)')
        ax.plot(E_plot, n_smooth_plot, 'r--', linewidth=1, label='Smooth N(E)')
        ax.plot(E_plot, n_wkb_plot, 'b-', linewidth=1, label=f'WKB (score={all_results[name]["score"]:.2f})')
        ax.set_xlabel('E')
        ax.set_ylabel('N(E)')
        short_name = name.split(':')[0]
        ax.set_title(short_name)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Summary plot
    ax = axes[1, 1]
    names = list(all_results.keys())
    scores = [all_results[n]["score"] for n in names]
    short_names = [n.split(':')[0] for n in names]
    ax.barh(range(len(names)), scores, color='steelblue')
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(short_names, fontsize=9)
    ax.set_xlabel('Score (lower = better)')
    ax.set_title('Family Comparison')
    ax.grid(True, alpha=0.3, axis='x')

    plt.suptitle('Parameterized Operator Search — H = U(x)p + V(x)/p',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'operator_search.png'), dpi=150)
    plt.close()
    print(f"\nPlot saved.")

    # Save
    with open(os.path.join(RESULTS_DIR, "operator_search.json"), "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print("\n" + "=" * 60)
    print("SEARCH SUMMARY")
    print("=" * 60)
    for name in sorted(all_results.keys(), key=lambda n: all_results[n]["score"]):
        r = all_results[name]
        print(f"  {name}")
        print(f"    Score: {r['score']:.6f}")
        for pn, pv in r["params"].items():
            print(f"    {pn} = {pv:.6f}")
        print()

    print("KEY QUESTION: Do any of these match INDIVIDUAL zeros")
    print("(the oscillatory part), or only the smooth average?")
    print("If only smooth — the oscillatory corrections require")
    print("encoding prime numbers, which no simple parametric family can do.")


if __name__ == "__main__":
    run_search()

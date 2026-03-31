"""
Deep heat kernel analysis of zeta zeros.

The heat trace Tr(e^{-tH}) = sum_n e^{-t*t_n} encodes ALL spectral information.
Its short-time asymptotics reveal the geometry of the underlying space:

  Tr(e^{-tH}) ~ sum_k a_k * t^((k-d)/2)  as t -> 0+

where d is the dimension and a_k are spectral invariants (Seeley-DeWitt coefficients).

For the zeta zeros, we found d ≈ 2.17 in the v2 analysis.
This script does a more careful analysis with the 2M dataset.
"""

import numpy as np
from scipy.optimize import curve_fit
import json
import os

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)


def load_zeros(source="lmfdb", count=None):
    if source == "2M":
        path = os.path.join(DATA_DIR, "zeros_odlyzko_2M")
    elif source == "lmfdb":
        path = os.path.join(DATA_DIR, "zeros_100k.txt")
    else:
        path = os.path.join(DATA_DIR, "zeros_odlyzko_100k")

    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if count and i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def heat_trace(zeros, t_values):
    """Compute sum_n exp(-t * t_n) for each t."""
    traces = []
    for t in t_values:
        # For very small t, many terms contribute
        # For large t, only the first few zeros matter
        tr = np.sum(np.exp(-t * zeros))
        traces.append(tr)
    return np.array(traces)


def analyze_short_time(zeros, label=""):
    """
    Analyze short-time heat kernel asymptotics.
    Tr(e^{-tH}) ~ a_0 * t^(-d/2) + a_1 * t^(-(d-2)/2) + ...

    We fit multiple models to determine d and the first few coefficients.
    """
    print(f"\n{'='*60}")
    print(f"HEAT KERNEL ANALYSIS — {label}")
    print(f"{'='*60}")
    print(f"Using {len(zeros)} zeros")

    # Compute heat trace at many t values
    t_values = np.logspace(-4, 1, 200)
    traces = heat_trace(zeros, t_values)

    # 1. Simple power law fit: Tr ~ A * t^alpha
    mask = t_values < 0.01  # use small t only
    log_t = np.log(t_values[mask])
    log_tr = np.log(traces[mask])

    coeffs = np.polyfit(log_t, log_tr, 1)
    alpha = coeffs[0]
    A = np.exp(coeffs[1])
    d_simple = -2 * alpha

    print(f"\n  Simple power law fit (t < 0.01):")
    print(f"    Tr(e^{{-tH}}) ~ {A:.4f} * t^({alpha:.4f})")
    print(f"    Effective dimension d = {d_simple:.4f}")

    # 2. Two-term fit: Tr ~ A*t^alpha + B*t^beta
    mask2 = t_values < 0.1
    t_fit = t_values[mask2]
    tr_fit = traces[mask2]

    def two_term(t, A, alpha, B, beta):
        return A * t**alpha + B * t**beta

    try:
        popt, pcov = curve_fit(two_term, t_fit, tr_fit,
                               p0=[A, alpha, 1.0, alpha + 0.5],
                               maxfev=10000)
        A2, alpha2, B2, beta2 = popt
        d_two = -2 * alpha2
        print(f"\n  Two-term fit (t < 0.1):")
        print(f"    Tr ~ {A2:.4f} * t^({alpha2:.4f}) + {B2:.4f} * t^({beta2:.4f})")
        print(f"    Leading dimension d = {d_two:.4f}")
        print(f"    Sub-leading exponent: {beta2:.4f}")
    except Exception as e:
        print(f"  Two-term fit failed: {e}")
        d_two = d_simple

    # 3. Compare with known models
    print(f"\n  Reference dimensions:")
    print(f"    d=1 (quantum mechanics on interval): alpha = -0.5")
    print(f"    d=2 (quantum mechanics on surface):  alpha = -1.0")
    print(f"    d=3 (3D quantum mechanics):          alpha = -1.5")
    print(f"    Our measured alpha:                   {alpha:.4f}")
    print(f"    Our measured d:                       {d_simple:.4f}")

    # 4. Weyl law check
    # For a d-dimensional operator, the counting function N(E) ~ C * E^(d/2)
    # For zeta zeros: N(t) ~ (t/2pi) * log(t/2pi) ~ t*log(t)/2pi
    # This grows faster than any power law, suggesting d is not a fixed integer
    # but the effective dimension changes with scale

    print(f"\n  Weyl law analysis:")
    print(f"    Zeta zero counting: N(t) ~ (t/2pi)*log(t/2pi)")
    print(f"    This is t*log(t) growth — FASTER than any t^(d/2)")
    print(f"    => The 'dimension' is scale-dependent (logarithmic correction)")
    print(f"    => The underlying space is not a manifold of fixed dimension")
    print(f"    => Suggests fractal or log-corrected geometry")

    # 5. Scale-dependent dimension
    print(f"\n  Scale-dependent effective dimension:")
    t_pairs = [(0.0001, 0.001), (0.001, 0.01), (0.01, 0.1), (0.1, 1.0)]
    for t1, t2 in t_pairs:
        tr1 = np.sum(np.exp(-t1 * zeros))
        tr2 = np.sum(np.exp(-t2 * zeros))
        local_alpha = np.log(tr2 / tr1) / np.log(t2 / t1)
        local_d = -2 * local_alpha
        print(f"    t in [{t1}, {t2}]: alpha={local_alpha:.4f}, d_eff={local_d:.4f}")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Log-log heat trace
    ax = axes[0, 0]
    ax.loglog(t_values, traces, 'b-', linewidth=1.5, label='Heat trace')
    t_line = t_values[t_values < 0.1]
    ax.loglog(t_line, A * t_line**alpha, 'r--', label=f'Fit: t^{{{alpha:.3f}}}')
    ax.set_xlabel('t')
    ax.set_ylabel('Tr(e^{-tH})')
    ax.set_title('Heat Kernel (log-log)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Local effective dimension
    ax = axes[0, 1]
    t_mid = []
    d_local = []
    for i in range(len(t_values) - 1):
        if traces[i] > 0 and traces[i+1] > 0:
            local_a = np.log(traces[i+1] / traces[i]) / np.log(t_values[i+1] / t_values[i])
            t_mid.append(np.sqrt(t_values[i] * t_values[i+1]))
            d_local.append(-2 * local_a)

    ax.semilogx(t_mid, d_local, 'b-', linewidth=1)
    ax.axhline(y=1, color='gray', linestyle=':', alpha=0.5, label='d=1')
    ax.axhline(y=2, color='gray', linestyle='--', alpha=0.5, label='d=2')
    ax.set_xlabel('t')
    ax.set_ylabel('Effective dimension d(t)')
    ax.set_title('Scale-Dependent Effective Dimension')
    ax.legend()
    ax.set_ylim(0, 4)
    ax.grid(True, alpha=0.3)

    # Heat trace * t^(d/2) — should approach a constant as t -> 0
    ax = axes[1, 0]
    normalized = traces * t_values**(d_simple/2)
    ax.semilogx(t_values, normalized, 'b-', linewidth=1.5)
    ax.set_xlabel('t')
    ax.set_ylabel(f'Tr(e^{{-tH}}) · t^{{{d_simple/2:.2f}}}')
    ax.set_title('Normalized Heat Trace (should plateau)')
    ax.grid(True, alpha=0.3)

    # Comparison: smooth counting function vs Weyl law
    ax = axes[1, 1]
    E_vals = np.linspace(10, 1000, 500)
    N_smooth = (E_vals / (2*np.pi)) * np.log(E_vals / (2*np.pi)) - E_vals / (2*np.pi) + 7/8
    N_weyl_d1 = 0.15 * E_vals  # d=1 Weyl law (linear)
    N_weyl_d2 = 0.005 * E_vals**1  # d=2 Weyl law

    ax.plot(E_vals, N_smooth, 'b-', linewidth=2, label='N(E) smooth')
    ax.plot(E_vals, N_weyl_d1, 'r--', alpha=0.5, label='Weyl d=1 (linear)')
    ax.set_xlabel('E')
    ax.set_ylabel('N(E)')
    ax.set_title('Zero Counting vs Weyl Law')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Heat Kernel Analysis — {label}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, f'heat_kernel_{label.lower().replace(" ", "_")}.png'), dpi=150)
    plt.close()
    print(f"\n  Plot saved.")

    return {
        "d_simple": float(d_simple),
        "alpha": float(alpha),
        "A": float(A),
        "scale_dependent_d": [(t1, t2, float(-2 * np.log(np.sum(np.exp(-t2*zeros)) / np.sum(np.exp(-t1*zeros))) / np.log(t2/t1))) for t1, t2 in t_pairs],
    }


def main():
    # Run on 100K
    zeros_100k = load_zeros("lmfdb", count=100000)
    r1 = analyze_short_time(zeros_100k, "100K LMFDB")

    # Run on 2M
    zeros_2M = load_zeros("2M")
    r2 = analyze_short_time(zeros_2M, "2M Odlyzko")

    # Compare
    print(f"\n{'='*60}")
    print("COMPARISON: 100K vs 2M")
    print(f"{'='*60}")
    print(f"  100K effective dimension: {r1['d_simple']:.4f}")
    print(f"  2M effective dimension:   {r2['d_simple']:.4f}")
    print(f"  Difference: {abs(r1['d_simple'] - r2['d_simple']):.4f}")

    print(f"\n  INTERPRETATION:")
    print(f"  The effective dimension is NOT a fixed integer.")
    print(f"  N(t) ~ t*log(t) means the Weyl law has a logarithmic correction,")
    print(f"  which is characteristic of operators on spaces with fractal")
    print(f"  or non-commutative geometry (cf. Connes' approach).")
    print(f"  This is a CONSTRAINT on the Hilbert-Polya operator:")
    print(f"  it cannot live on an ordinary d-dimensional manifold.")

    results = {"100k": r1, "2M": r2}
    with open(os.path.join(RESULTS_DIR, "heat_kernel_analysis.json"), "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fit the interior margin (λ_min) scaling law vs conductor c.

Tests three models:
1. Power law: λ_min ~ c^(-α)
2. Exponential: λ_min ~ exp(-β * c)
3. Stretched exponential: λ_min ~ exp(-β * c^γ)

Also checks the coincidence:
- c ≤ 100 ⟺ a ≤ ½log(100) ≈ 2.30
- Largest verified zero height T ≈ 3e12 gives a_max ≈ log(T)/(4π) ≈ 2.29
"""

import json
import sys
import math
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import xlogy
from pathlib import Path

def fit_scaling(results_file, output_prefix="scaling_fit"):
    """Fit scaling laws to the data."""

    # Load results
    with open(results_file, 'r') as f:
        results = json.load(f)

    # Extract data, converting log10 to natural log
    c_vals = []
    lmin_vals = []
    log_lmin_vals = []

    for r in results:
        if isinstance(r, dict):
            c = r['c']
            lmin = r['lambda_min']
        else:
            continue

        if lmin <= 0:
            continue

        c_vals.append(c)
        lmin_vals.append(lmin)
        log_lmin_vals.append(math.log(lmin))

    c_vals = np.array(c_vals)
    lmin_vals = np.array(lmin_vals)
    log_lmin_vals = np.array(log_lmin_vals)

    print(f"=== SCALING LAW FITTING ===")
    print(f"Data points: {len(c_vals)}")
    print(f"c range: {c_vals.min()} - {c_vals.max()}")
    print(f"λ_min range: {lmin_vals.min():.6e} - {lmin_vals.max():.6e}")
    print()

    fits = {}

    # Model 1: Power law λ_min = A * c^(-α)
    print("--- MODEL 1: Power Law λ_min = A * c^(-α) ---")
    try:
        def power_law(c, A, alpha):
            return A * np.power(c, -alpha)

        # Initial guess
        popt, pcov = curve_fit(power_law, c_vals, lmin_vals, p0=[1, 5], maxfev=5000)
        A, alpha = popt

        # Compute R²
        y_pred = power_law(c_vals, A, alpha)
        ss_res = np.sum((lmin_vals - y_pred) ** 2)
        ss_tot = np.sum((lmin_vals - lmin_vals.mean()) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        print(f"  A = {A:.6e}")
        print(f"  α = {alpha:.4f}")
        print(f"  R² = {r2:.6f}")
        print(f"  Fit: λ_min = {A:.6e} * c^({-alpha:.4f})")
        print()

        fits['power_law'] = {
            'A': A,
            'alpha': alpha,
            'r2': r2,
            'model': 'lambda_min = A * c^(-alpha)'
        }
    except Exception as e:
        print(f"  ERROR: {e}\n")

    # Model 2: Exponential decay λ_min = A * exp(-β * c)
    print("--- MODEL 2: Exponential λ_min = A * exp(-β * c) ---")
    try:
        def exp_decay(c, A, beta):
            return A * np.exp(-beta * c)

        # Initial guess
        popt, pcov = curve_fit(exp_decay, c_vals, lmin_vals, p0=[0.1, 0.3], maxfev=5000)
        A, beta = popt

        # Compute R²
        y_pred = exp_decay(c_vals, A, beta)
        ss_res = np.sum((lmin_vals - y_pred) ** 2)
        ss_tot = np.sum((lmin_vals - lmin_vals.mean()) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        print(f"  A = {A:.6e}")
        print(f"  β = {beta:.6f}")
        print(f"  R² = {r2:.6f}")
        print(f"  Fit: λ_min = {A:.6e} * exp({-beta:.6f} * c)")
        print()

        fits['exponential'] = {
            'A': A,
            'beta': beta,
            'r2': r2,
            'model': 'lambda_min = A * exp(-beta * c)'
        }
    except Exception as e:
        print(f"  ERROR: {e}\n")

    # Model 3: Stretched exponential λ_min = A * exp(-β * c^γ)
    print("--- MODEL 3: Stretched Exponential λ_min = A * exp(-β * c^γ) ---")
    try:
        def stretched_exp(c, A, beta, gamma):
            return A * np.exp(-beta * np.power(c, gamma))

        # Initial guess
        popt, pcov = curve_fit(stretched_exp, c_vals, lmin_vals, p0=[0.1, 0.1, 1.5], maxfev=10000)
        A, beta, gamma = popt

        # Compute R²
        y_pred = stretched_exp(c_vals, A, beta, gamma)
        ss_res = np.sum((lmin_vals - y_pred) ** 2)
        ss_tot = np.sum((lmin_vals - lmin_vals.mean()) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        print(f"  A = {A:.6e}")
        print(f"  β = {beta:.6f}")
        print(f"  γ = {gamma:.4f}")
        print(f"  R² = {r2:.6f}")
        print(f"  Fit: λ_min = {A:.6e} * exp({-beta:.6f} * c^{gamma:.4f})")
        print()

        fits['stretched_exp'] = {
            'A': A,
            'beta': beta,
            'gamma': gamma,
            'r2': r2,
            'model': 'lambda_min = A * exp(-beta * c^gamma)'
        }
    except Exception as e:
        print(f"  ERROR: {e}\n")

    # Find best fit
    best_model = max(fits.items(), key=lambda x: x[1]['r2'])
    print(f"=== BEST FIT: {best_model[0]} (R² = {best_model[1]['r2']:.6f}) ===")
    print()

    # Compute boundary coincidence check
    print("=== BOUNDARY COINCIDENCE CHECK ===")
    c_edge = 100
    a_edge = 0.5 * math.log(c_edge)
    print(f"Theorem edge: c ≤ {c_edge} ⟺ a ≤ ½log({c_edge}) = {a_edge:.4f}")

    T_max = 3e12  # Platt's largest verified zero height
    a_max_platt = math.log(T_max) / (4 * math.pi)
    print(f"Platt maximum: largest verified zero at height T ≈ {T_max:.2e}")
    print(f"  a_max = log(T)/(4π) = log({T_max:.2e})/(4π) = {a_max_platt:.4f}")

    diff = abs(a_edge - a_max_platt)
    rel_diff = diff / a_edge

    print()
    print(f"Comparison:")
    print(f"  a_edge = {a_edge:.6f}")
    print(f"  a_max  = {a_max_platt:.6f}")
    print(f"  Δa = {diff:.6f}")
    print(f"  Relative difference: {rel_diff * 100:.3f}%")

    if rel_diff < 0.01:  # Less than 1% difference
        print(f"✓ COINCIDENCE DETECTED: These match to within {rel_diff*100:.2f}% — structural!")
    elif rel_diff < 0.05:  # Less than 5%
        print(f"~ Possible coincidence: {rel_diff*100:.2f}% difference")
    else:
        print(f"✗ No coincidence detected: {rel_diff*100:.2f}% difference")

    # Save results
    output_file = Path(output_prefix + "_results.json")
    with open(output_file, 'w') as f:
        json.dump({
            "data": {
                "c": c_vals.tolist(),
                "lambda_min": lmin_vals.tolist(),
                "log_lambda_min": log_lmin_vals.tolist()
            },
            "fits": fits,
            "best_fit": {
                "model": best_model[0],
                "r2": best_model[1]['r2'],
                "parameters": {k: v for k, v in best_model[1].items() if k not in ['r2', 'model']}
            },
            "boundary_check": {
                "c_edge": c_edge,
                "a_edge": a_edge,
                "T_max": T_max,
                "a_max_platt": a_max_platt,
                "delta_a": diff,
                "relative_diff": rel_diff,
                "coincident": rel_diff < 0.01
            }
        }, f, indent=2)

    print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: fit_scaling_law.py <results.json> [output_prefix]")
        sys.exit(1)

    results_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else "scaling_fit"

    fit_scaling(results_file, output_prefix)

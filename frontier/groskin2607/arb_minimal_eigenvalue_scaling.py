#!/usr/bin/env python3
"""
Extract the minimum eigenvalue via Arb high precision and fit scaling laws.

This is similar to margin_law.py but with focus on understanding the 0.44 exponent
through the lens of matrix structure and perturbation theory.

Key hypothesis to test:
  λ_min(c,N) ~ K · c^(-α(N)·N) · f(c,N)

where α(N) might be constant (α ≈ 0.44) or might depend on N.
"""

import sys
import argparse
import numpy as np
import json
import math
from pathlib import Path

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')

from flint import arb, ctx
import mpmath as mp


def extract_min_eigenvalue_arb(c, N, prec):
    """
    Extract λ_min by computing eigenvalues via mpmath, using Arb matrix.
    """
    from arb_ldlt_certify import build_arb_tau

    ctx.prec = prec
    A_arb, DIM = build_arb_tau(c, N, prec)

    # Convert to mpmath for eigenvalue computation
    mp_dps = int(prec * 0.29)
    mp.mp.dps = mp_dps

    A_mp = mp.matrix(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            A_mp[i, j] = mp.mpf(A_arb[i, j].str(int(prec*0.28), radius=False))

    evals = sorted([float(e) for e in mp.eigsy(A_mp, eigvals_only=True)])
    return evals[0]  # minimum eigenvalue


def scan_over_c_at_fixed_N(N, c_values, prec=1000):
    """
    Scan λ_min over varying c at fixed N.

    Returns:
      (c_array, lambda_array, log_c_array, log_lambda_array)
    """
    print(f"\nScanning c at fixed N={N}...")
    lambda_vals = []

    for c in c_values:
        try:
            lam = extract_min_eigenvalue_arb(c, N, prec)
            lambda_vals.append(lam)
            print(f"  c={c:5d}: λ_min = {lam:+.6e}  log₁₀(|λ|) = {np.log10(abs(lam)+1e-20):+.2f}")
        except Exception as e:
            print(f"  c={c:5d}: FAILED ({e})")
            lambda_vals.append(None)

    c_array = np.array(c_values)
    lambda_array = np.array(lambda_vals)
    log_c_array = np.log10(c_array)
    log_lambda_array = np.array([np.log10(abs(l) + 1e-20) if l is not None else None for l in lambda_vals])

    return c_array, lambda_array, log_c_array, log_lambda_array


def fit_power_law(log_c, log_lambda):
    """
    Fit log(λ) = a·log(c) + b, returning the exponent a.
    """
    # Remove None values
    valid = [(lc, ll) for lc, ll in zip(log_c, log_lambda) if ll is not None]
    if len(valid) < 2:
        return None
    log_c_valid, log_lambda_valid = zip(*valid)
    log_c_valid = np.array(log_c_valid)
    log_lambda_valid = np.array(log_lambda_valid)

    coeffs = np.polyfit(log_c_valid, log_lambda_valid, 1)
    return coeffs[0]  # slope


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--N-values', type=str, default='4,6,8', help='comma-separated N values')
    p.add_argument('--c-values', type=str, default='13,20,30,50', help='comma-separated c values')
    p.add_argument('--prec', type=int, default=1000, help='Arb precision in bits')
    args = p.parse_args()

    N_vals = [int(x) for x in args.N_values.split(',')]
    c_vals = [int(x) for x in args.c_values.split(',')]

    print("="*100)
    print("ARB MINIMAL EIGENVALUE SCALING ANALYSIS")
    print("="*100)

    results = {}

    # Scan over c at each fixed N
    for N in N_vals:
        print(f"\n{'='*100}")
        print(f"N = {N}")
        print(f"{'='*100}")

        c_array, lambda_array, log_c_array, log_lambda_array = scan_over_c_at_fixed_N(
            N, c_vals, prec=args.prec
        )

        # Fit power law
        exponent = fit_power_law(log_c_array, log_lambda_array)
        print(f"\nPower-law fit (over c):")
        if exponent is not None:
            print(f"  log₁₀(λ) = {exponent:.4f}·log₁₀(c) + b")
            print(f"  → λ ~ c^{exponent:.4f}")
            print(f"  Exponent/N = {exponent / N:.4f}  (should be -0.44 if hypothesis holds)")

        results[f'N={N}'] = {
            'c_values': c_vals,
            'lambda_values': [float(x) if x is not None else None for x in lambda_array],
            'exponent': exponent,
        }

    # Scan over N at each fixed c
    print(f"\n{'='*100}")
    print("Scanning N at fixed c")
    print(f"{'='*100}")

    for c in c_vals:
        print(f"\nAt c={c}:")
        lambda_at_c = []

        for N in N_vals:
            try:
                lam = extract_min_eigenvalue_arb(c, N, args.prec)
                lambda_at_c.append(lam)
                print(f"  N={N}: λ_min = {lam:+.6e}  log₁₀(|λ|) = {np.log10(abs(lam)+1e-20):+.2f}")
            except Exception as e:
                print(f"  N={N}: FAILED")
                lambda_at_c.append(None)

        # Fit quadratic in N (as margin_law.py did)
        valid_N = [N_vals[i] for i in range(len(N_vals)) if lambda_at_c[i] is not None]
        valid_lam = [np.log10(abs(l) + 1e-20) for l in lambda_at_c if l is not None]

        if len(valid_N) >= 2:
            quad_fit = np.polyfit(valid_N, valid_lam, 2)
            print(f"\n  Quadratic fit in N:")
            print(f"    log₁₀(λ) = {quad_fit[0]:.4f}·N² + {quad_fit[1]:.4f}·N + {quad_fit[2]:.2f}")
            print(f"    Leading term suggests scaling ~ c^({quad_fit[1]:.4f}/log(c)·N) ≈ c^({quad_fit[1]/(np.log10(c)):.4f}·N)")

        results[f'c={c}'] = {
            'N_values': N_vals,
            'lambda_values': [float(x) if x is not None else None for x in lambda_at_c],
        }

    # Save
    output_json = 'arb_minimal_eigenvalue_scaling.json'
    with open(output_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n\nSaved results to {output_json}")

    print("\n" + "="*100)
    print("SUMMARY")
    print("="*100)
    print("\nTo verify the c^(-0.44N) scaling:")
    print("  - Exponent/N should be close to -0.44")
    print("  - If all exponents are around 0.4-0.5, the universality class might be")
    print("    related to Toeplitz-plus-low-rank matrices")


if __name__ == '__main__':
    main()

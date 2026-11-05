#!/usr/bin/env python3
"""
Zero-fitting optimization for CCM approximation bound.

Given the computed ground-state eigenvector xi_0, find the perturbation delta_xi
that best improves the fit to true Riemann zeta zeros, subject to remaining an
eigenvector (or at least, improving zero fit).

This is a Fourier-coefficient least-squares approach, not a full operator
reconstruction.
"""

import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')

from connes_decisive_v3 import build_tau, xi_hat_factory, find_zeros, TRUE_ZEROS, L, TWO_PI
import numpy as np
from scipy.optimize import minimize, differential_evolution
import time

def zero_fit_error(xi_flat, idx, target_zeros, xh_factory=xi_hat_factory, verbose=False):
    """
    Compute sum of squared distances from computed zeros to target zeros.

    Args:
        xi_flat: flattened coefficient vector
        idx: index array (pole positions)
        target_zeros: list of target zero locations
        xh_factory: function to build Fourier function from coefficients
        verbose: if True, print debug info

    Returns:
        error: sum of squared distances
    """
    xi = np.array(xi_flat, dtype=float)

    # Normalize to match original xi_0 (avoid trivial scaling)
    s = np.sum(xi) / np.sqrt(L)
    if abs(s) > 1e-30:
        xi = xi / s

    try:
        xh = xh_factory(xi, idx)
        zeros = find_zeros(xh, 10, 60, 2000)
    except:
        return 1e10  # Return large error on failure

    if not zeros:
        return 1e10

    # Compute sum of squared distances to nearest target zero
    error = 0.0
    matched_count = 0
    for tz in target_zeros:
        nearest = min(zeros, key=lambda z: abs(z - tz))
        dist = abs(nearest - tz)
        error += dist ** 2
        matched_count += 1
        if verbose:
            print(f"  Target {tz:.6f}: nearest {nearest:.6f}, dist {dist:.2e}")

    if verbose:
        print(f"  Total error: {error:.4e} (matched {matched_count} targets)")

    return error


def measure_zero_fit(xi, idx, target_zeros, label=""):
    """Helper: compute and report zero-fit error for a given xi."""
    error = zero_fit_error(xi.flatten(), idx, target_zeros)
    return error


def run_optimization(N, target_zeros=None):
    """
    Build tau at size N, extract ground-state eigenvector,
    optimize for zero-fit improvement, and measure perturbation.
    """
    if target_zeros is None:
        target_zeros = TRUE_ZEROS[:6]  # Default: first 6 true zeros

    print(f"\n{'='*70}")
    print(f"N={N}: Building tau and extracting eigenvector")
    print(f"{'='*70}")

    t0 = time.time()
    R = TWO_PI * N / L + 80.0
    G = max(60000, 3000 * N)

    # Build operator
    tau, idx = build_tau(N, R, G)
    print(f"Tau built in {time.time()-t0:.1f}s, shape {tau.shape}")

    # Eigenvector extraction
    t0 = time.time()
    ev, evec = np.linalg.eigh(tau)
    k0 = int(np.argmin(np.abs(ev)))
    v = evec[:, k0]
    D = 2 * N + 1

    # Extract even part
    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])
    s = even.sum() / np.sqrt(L)
    if abs(s) > 1e-30:
        even = even / s

    xi_0 = even.copy()
    print(f"Eigenvector extracted in {time.time()-t0:.1f}s")
    print(f"  eig[0]={ev[0]:.4e}, eig[k0]={ev[k0]:.4e}")
    print(f"  xi_0 shape={xi_0.shape}, norm={np.linalg.norm(xi_0):.4f}")

    # Compute initial zero-fit error
    print(f"\nComputing initial zero fit (target zeros: {len(target_zeros)})")
    error_0 = zero_fit_error(xi_0.flatten(), idx, target_zeros, verbose=True)
    print(f"Initial error: {error_0:.4e}")

    # Optimize for better zero fit
    print(f"\nOptimizing zero fit...")
    t_opt = time.time()

    # Use differential evolution for global search (more robust than local minimize)
    bounds = [(x - 0.5, x + 0.5) for x in xi_0]  # Search near xi_0
    result = differential_evolution(
        lambda xf: zero_fit_error(xf, idx, target_zeros),
        bounds,
        seed=42,
        maxiter=300,
        popsize=30,
        atol=1e-6,
        tol=1e-6,
        workers=1,
        updating='deferred'
    )

    opt_time = time.time() - t_opt
    print(f"Optimization completed in {opt_time:.1f}s, {result.nfev} function evals")
    print(f"Result: success={result.success}, nfev={result.nfev}, fun={result.fun:.4e}")

    xi_opt = result.x.copy()

    # Renormalize
    s_opt = np.sum(xi_opt) / np.sqrt(L)
    if abs(s_opt) > 1e-30:
        xi_opt = xi_opt / s_opt

    error_opt = zero_fit_error(xi_opt.flatten(), idx, target_zeros, verbose=True)
    print(f"Optimized error: {error_opt:.4e}")

    # Measure perturbation
    delta_xi = xi_opt - xi_0
    epsilon = np.linalg.norm(delta_xi)

    print(f"\n{'='*70}")
    print(f"RESULTS for N={N}:")
    print(f"{'='*70}")
    print(f"Initial zero-fit error:    {error_0:.4e}")
    print(f"Optimized zero-fit error:  {error_opt:.4e}")
    print(f"Improvement:               {error_0 - error_opt:.4e} ({100*(error_0-error_opt)/max(error_0,1e-10):.1f}%)")
    print(f"Eigenvector perturbation:  ε = {epsilon:.4e}")
    print(f"Relative perturbation:     ε/||xi_0|| = {epsilon/np.linalg.norm(xi_0):.4e}")

    return {
        'N': N,
        'error_initial': error_0,
        'error_optimized': error_opt,
        'improvement': error_0 - error_opt,
        'epsilon': epsilon,
        'rel_epsilon': epsilon / np.linalg.norm(xi_0),
        'xi_0': xi_0,
        'xi_opt': xi_opt,
        'delta_xi': delta_xi,
        'idx': idx,
        'target_zeros': target_zeros,
    }


if __name__ == "__main__":
    # Run at N=32 (main result)
    result_n32 = run_optimization(N=32, target_zeros=TRUE_ZEROS[:6])

    # Save results
    import json
    output = {
        'N': result_n32['N'],
        'target_zeros': result_n32['target_zeros'],
        'error_initial': float(result_n32['error_initial']),
        'error_optimized': float(result_n32['error_optimized']),
        'improvement': float(result_n32['improvement']),
        'epsilon': float(result_n32['epsilon']),
        'rel_epsilon': float(result_n32['rel_epsilon']),
    }

    with open('/Users/scones/projects/riemann/results/zero_fit_optimization_n32.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to results/zero_fit_optimization_n32.json")


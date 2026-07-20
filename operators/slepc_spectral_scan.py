#!/usr/bin/env python3
"""
Large-scale spectral scan using scipy.sparse.linalg.eigsh (ARPACK).

This script rebuilds the Connes-Consani operator at increasing matrix sizes
and scans for convergence of eigenvalues to known Riemann zeta zeros.

The operator construction is identical to connes_v5_spectrum.py:
- Weil matrix τ = arch + B - P (archimedean + boundary + prime terms)
- Extract minimum eigenvector ξ
- Find zeros of P(s) = Σ_j ξ_j/(j-s)
- Physical spectrum: γ = 2π·s/L
"""

import os
import sys
import json
import time
import math
import numpy as np
from scipy import sparse
from scipy.sparse import linalg as splinalg
from scipy.integrate import quad
from scipy.optimize import brentq

import mpmath

DATA_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

mpmath.mp.dps = 80

LAMBDA_SQ = 14
L         = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))
PRIMES    = [2, 3, 5, 7, 11, 13]


def load_zeros(count=100):
    """Load reference zeros from data file."""
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    try:
        with open(path) as f:
            for i, line in enumerate(f):
                if i >= count:
                    break
                line = line.strip()
                if line:
                    zeros.append(float(line))
    except FileNotFoundError:
        print(f"WARNING: Could not load zeros from {path}")
    return np.array(zeros)


def rho_scalar(x):
    """ρ(x) = e^{x/2} / (2sinh(x))"""
    if x < 1e-10:
        return 1.0 / (2 * x + 1e-300)
    return np.exp(x / 2) / (np.exp(x) - np.exp(-x))


def alpha_L(n, L):
    """α_L(n) = (1/π) ∫₀ᴸ sin(2πnx/L) ρ(x) dx"""
    if n == 0:
        return 0.0
    a = 2 * np.pi * n / L

    def integrand(x):
        if x < 1e-12:
            return a / 2.0
        return np.sin(a * x) * rho_scalar(x)

    result, _ = quad(integrand, 0, L, limit=300, epsabs=1e-12, epsrel=1e-12,
                     points=[1e-8, 0.01, 0.1, L / 2])
    return result / np.pi


def W_arch_diag(n, L):
    """W_R(V_n, V_n) = Re[ψ(1/4 + iπn/L)] + log(2)"""
    if n == 0:
        return float(-mpmath.euler - mpmath.log(2))
    s = float(mpmath.pi) * n / L
    return float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j * s))) + float(mpmath.log(2))


def W_boundary(n, m, L):
    """Boundary term in Weil matrix."""
    pi   = float(mpmath.pi)
    L_mp = mpmath.mpf(str(L))
    num  = 32 * L_mp * mpmath.sinh(L_mp / 4)**2 * (L_mp**2 - 16 * pi**2 * m * n)
    den  = (L_mp**2 + 16 * pi**2 * m**2) * (L_mp**2 + 16 * pi**2 * n**2)
    if abs(float(den)) < 1e-50:
        return 0.0
    return float(num / den)


def von_mangoldt(k, primes):
    """von Mangoldt function."""
    if k <= 1: return 0.0
    for p in primes:
        if p > k: break
        pk = p
        while pk <= k:
            if pk == k: return float(mpmath.log(p))
            pk *= p
    return 0.0


def W_prime(n, m, L, primes):
    """Prime term in Weil matrix."""
    pi    = float(mpmath.pi)
    k_max = int(float(mpmath.exp(mpmath.mpf(str(L))))) + 1
    total = 0.0
    for k in range(2, k_max + 1):
        lam_k = von_mangoldt(k, primes)
        if lam_k == 0: continue
        y      = float(mpmath.log(k))
        weight = lam_k * k**(-0.5)
        if n != m:
            q = (np.sin(2 * pi * m * y / L) - np.sin(2 * pi * n * y / L)) / (pi * (n - m))
        else:
            q = 2 * (1 - abs(y) / L) * np.cos(2 * pi * n * y / L)
        total += weight * q
    return total


def build_weil_matrix(N, L, primes, alpha_cache):
    """Build the Weil matrix."""
    dim = 2 * N + 1
    tau = np.zeros((dim, dim))

    for i in range(dim):
        n = i - N
        for j in range(i, dim):
            m = j - N

            # Archimedean
            if n == m:
                arch = W_arch_diag(n, L)
            else:
                a_m = alpha_cache[m]
                a_n = alpha_cache[n]
                arch = (a_m - a_n) / (n - m)

            # Boundary
            bnd = W_boundary(n, m, L)

            # Prime
            prm = W_prime(n, m, L, primes)

            val = arch + bnd - prm
            tau[i, j] = val
            tau[j, i] = val

    return tau


def find_spectrum_rational(xi_vec, N, s_min=-80, s_max=80):
    """Find zeros of P(s) = Σ_{j=-N}^{N} ξ_j / (j - s)"""
    js = np.arange(-N, N + 1, dtype=float)

    def P(s):
        diffs = js - s
        safe  = np.where(np.abs(diffs) < 1e-12, 1e-12, diffs)
        return float(np.sum(xi_vec / safe))

    zeros = []
    for j_int in range(int(s_min), int(s_max)):
        lo = j_int + 1e-8
        hi = j_int + 1 - 1e-8
        try:
            f_lo = P(lo)
            f_hi = P(hi)
            if f_lo * f_hi < 0:
                z = brentq(P, lo, hi, xtol=1e-12, rtol=1e-12)
                zeros.append(z)
        except Exception:
            pass

    return np.array(zeros)


def compute_spectrum_at_N(N, L, primes, ref_zeros, verbose=True):
    """
    Compute spectrum at a given matrix size N.

    Returns:
        dict with keys: N, time_alpha, time_matrix, time_eigen, n_match_tight,
                        s_zeros, gamma_pred, min_eigenvalue, errors
    """
    result = {"N": N}

    if verbose:
        print(f"\n{'='*60}")
        print(f"Computing at N = {N} (matrix dim: {2*N+1}×{2*N+1})")
        print(f"{'='*60}")

    # Compute α_L
    if verbose:
        print(f"Computing α_L(n) for n = -{N}..{N}...")
    t0 = time.time()
    alpha_cache = {}
    for n in range(-N, N + 1):
        alpha_cache[n] = alpha_L(n, L)
    time_alpha = time.time() - t0
    result["time_alpha"] = time_alpha
    if verbose:
        print(f"  Done in {time_alpha:.2f}s")

    # Build Weil matrix
    if verbose:
        print(f"Building {2*N+1}×{2*N+1} Weil matrix...")
    t0 = time.time()
    tau = build_weil_matrix(N, L, primes, alpha_cache)
    time_matrix = time.time() - t0
    result["time_matrix"] = time_matrix
    if verbose:
        print(f"  Built in {time_matrix:.2f}s")

    # Compute eigenvalues using scipy.sparse.linalg.eigsh (ARPACK)
    # We need the smallest eigenvalue and its eigenvector
    if verbose:
        print(f"Computing eigenvalues using ARPACK (scipy.sparse.linalg)...")
    t0 = time.time()

    # Convert to sparse for eigsh
    tau_sparse = sparse.csr_matrix(tau)

    # Use which='SA' for smallest algebraic values
    try:
        evals, evecs = splinalg.eigsh(tau_sparse, k=min(5, 2*N), which='SA',
                                      tol=1e-10, maxiter=1000)
    except Exception as e:
        if verbose:
            print(f"  ARPACK failed: {e}, trying dense fallback...")
        evals, evecs = np.linalg.eigh(tau)

    time_eigen = time.time() - t0
    result["time_eigen"] = time_eigen
    if verbose:
        print(f"  Done in {time_eigen:.2f}s")

    eps = float(evals[0])
    result["min_eigenvalue"] = eps
    xi = evecs[:, 0]

    if verbose:
        print(f"ε_N = {eps:.8e}")

    # Normalize and find spectrum
    eta_xi = np.sum(xi)
    if abs(eta_xi) > 1e-12:
        xi_norm = xi / eta_xi
    else:
        xi_norm = xi

    if verbose:
        print(f"Finding spectrum zeros...")
    s_zeros = find_spectrum_rational(xi_norm, N, s_min=0, s_max=N)
    gamma_pred = 2 * np.pi * s_zeros / L

    if verbose:
        print(f"Found {len(s_zeros)} zeros")

    # Match to reference zeros
    n_ref = min(50, len(ref_zeros))
    n_pred = min(n_ref, len(gamma_pred))
    tight_tol = 0.5
    n_match_tight = 0
    errors = []

    for i in range(n_pred):
        ref = ref_zeros[i]
        pred = gamma_pred[i]
        err = abs(pred - ref)
        errors.append(err)
        if err < tight_tol:
            n_match_tight += 1

    result["n_match_tight"] = n_match_tight
    result["n_compared"] = n_pred
    result["s_zeros"] = s_zeros[:n_pred].tolist() if len(s_zeros) > 0 else []
    result["gamma_pred"] = gamma_pred[:n_pred].tolist() if len(gamma_pred) > 0 else []
    result["errors"] = errors
    result["mean_error"] = float(np.mean(errors)) if len(errors) > 0 else None
    result["max_error"] = float(np.max(errors)) if len(errors) > 0 else None

    if verbose:
        print(f"Matched tight (<{tight_tol}): {n_match_tight}/{n_pred}")
        if len(errors) > 0:
            print(f"Mean error: {np.mean(errors):.6e}, Max error: {np.max(errors):.6e}")

    return result


def main():
    print("\n" + "="*70)
    print("SLEPC/ARPACK SPECTRAL SCAN — Convergence to Riemann Zeros")
    print("="*70)
    print(f"Solver: scipy.sparse.linalg.eigsh (ARPACK backend)")
    print(f"λ²={LAMBDA_SQ}, L={L:.6f}")

    # Load reference zeros
    ref_zeros = load_zeros(100)
    print(f"Loaded {len(ref_zeros)} reference Riemann zeros")

    # Test sizes: 500 (sanity check), then scale up
    test_sizes = [500, 1000, 2000, 5000]

    results = {}
    all_results = {
        "solver": "scipy.sparse.linalg.eigsh (ARPACK)",
        "lambda_sq": LAMBDA_SQ,
        "L": L,
        "scans": []
    }

    for N in test_sizes:
        try:
            result = compute_spectrum_at_N(N, L, PRIMES, ref_zeros, verbose=True)
            results[N] = result
            all_results["scans"].append(result)

            # Quick convergence check
            if len(all_results["scans"]) > 1:
                prev = all_results["scans"][-2]
                curr = result
                if curr["mean_error"] is not None and prev["mean_error"] is not None:
                    improvement = (prev["mean_error"] - curr["mean_error"]) / prev["mean_error"] * 100
                    print(f"\nError improvement vs previous: {improvement:.1f}%")

        except Exception as e:
            print(f"\nERROR at N={N}: {e}")
            import traceback
            traceback.print_exc()
            all_results["scans"].append({"N": N, "error": str(e)})

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for scan in all_results["scans"]:
        if "error" not in scan:
            print(f"N={scan['N']:5d}: {scan['n_match_tight']:2d}/{scan['n_compared']} matches, " +
                  f"mean_err={scan['mean_error']:.6e}, max_err={scan['max_error']:.6e}")
        else:
            print(f"N={scan['N']:5d}: FAILED - {scan['error']}")

    # Save results
    out_file = os.path.join(RESULTS_DIR, "slepc_spectral_scan.json")
    with open(out_file, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_file}")


if __name__ == "__main__":
    main()

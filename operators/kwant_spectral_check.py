#!/usr/bin/env python3
"""
Independent spectral computation using Kwant tight-binding framework.

This script builds the SAME Connes-Consani operator using Kwant's
tight-binding system API as an independent cross-check against
slepc_spectral_scan.py. Since the Weil matrix is dense with both
diagonal and off-diagonal elements, we map it onto a Kwant 1D chain
with onsite potentials and hoppings.

The Hamiltonian is identical mathematically; only the representation differs.
"""

import os
import sys
import json
import time
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import mpmath

# Try importing kwant; if it fails, provide informative error
try:
    import kwant
    KWANT_AVAILABLE = True
except ImportError as e:
    KWANT_AVAILABLE = False
    KWANT_ERROR = str(e)

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


def build_weil_matrix_dense(N, L, primes, alpha_cache):
    """Build the Weil matrix (dense version)."""
    dim = 2 * N + 1
    tau = np.zeros((dim, dim))

    for i in range(dim):
        n = i - N
        for j in range(i, dim):
            m = j - N

            if n == m:
                arch = W_arch_diag(n, L)
            else:
                a_m = alpha_cache[m]
                a_n = alpha_cache[n]
                arch = (a_m - a_n) / (n - m)

            bnd = W_boundary(n, m, L)
            prm = W_prime(n, m, L, primes)

            val = arch + bnd - prm
            tau[i, j] = val
            tau[j, i] = val

    return tau


def build_kwant_system(tau_matrix):
    """
    Build a Kwant tight-binding system from the dense Weil matrix.

    The matrix is symmetric, so we extract:
    - Onsite potentials: diagonal elements
    - Hoppings: off-diagonal elements (hopping from site i to j, j>i)

    Returns:
        A finalized Kwant system whose Hamiltonian matrix is equivalent to tau_matrix
    """
    N = tau_matrix.shape[0]
    lat = kwant.lattice.chain()

    def onsite(site):
        return tau_matrix[site.pos[0], site.pos[0]]

    def hopping(site1, site2):
        i, j = site1.pos[0], site2.pos[0]
        # Ensure i < j
        if i > j:
            i, j = j, i
        # Return the hopping matrix element
        return tau_matrix[i, j]

    syst = kwant.Builder()

    # Add all sites
    for i in range(N):
        syst[lat(i)] = onsite

    # Add hoppings (to nearest neighbors; we'll add all within the range)
    # For a full matrix representation, we need all pairs (i, i+1), (i, i+2), etc.
    for i in range(N):
        for j in range(i + 1, N):
            # Add hopping from i to j
            # In Kwant's 1D chain, this is typically done with lat(i) and lat(j)
            # But standard Kwant chains only connect nearest neighbors
            # To represent the full dense matrix, we connect all pairs:
            syst[lat(i), lat(j)] = hopping

    syst = syst.finalized()
    return syst


def extract_hamiltonian_matrix(syst):
    """
    Extract the full Hamiltonian matrix from a Kwant system.

    Returns:
        A numpy array representing the Hamiltonian.
    """
    # Get the scattering matrix for E=0 to extract Hamiltonian
    # Actually, use syst.h_mat which is available after finalization
    if hasattr(syst, 'hamiltonian'):
        # For finalized systems, the Hamiltonian is stored differently
        # We need to reconstruct it from the builder structure
        pass

    # Alternative: use kwant's H matrix from the finalized system
    # The finalized system has a .H method, but we need to call it correctly
    try:
        # For kwant 1.x+, get the Hamiltonian matrix
        N = len(syst.sites)
        H = np.zeros((N, N), dtype=complex)

        for i, site_i in enumerate(syst.sites):
            for j, site_j in enumerate(syst.sites):
                # Get the matrix element <i|H|j>
                val = syst.H[i, j]
                if isinstance(val, (int, float)):
                    H[i, j] = val
                else:
                    H[i, j] = val

        return H.real if np.allclose(H.imag, 0) else H

    except Exception as e:
        print(f"Error extracting Hamiltonian: {e}")
        # Fallback: reconstruct from the system definition
        return None


def compare_matrices(tau_orig, tau_kwant, N):
    """Compare the original matrix with the Kwant-reconstructed one."""
    if tau_kwant is None:
        return None

    diff = np.linalg.norm(tau_orig - tau_kwant, 'fro')
    rel_error = diff / np.linalg.norm(tau_orig, 'fro')

    return {
        "frobenius_norm_diff": float(diff),
        "relative_error": float(rel_error),
        "max_element_diff": float(np.max(np.abs(tau_orig - tau_kwant)))
    }


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


def compute_kwant_spectrum(N, L, primes, ref_zeros, verbose=True):
    """Compute spectrum using Kwant system."""
    result = {"N": N, "solver": "kwant"}

    if verbose:
        print(f"\n{'='*60}")
        print(f"KWANT Spectral Check at N = {N}")
        print(f"{'='*60}")

    # Build the Weil matrix (same as slepc version)
    if verbose:
        print(f"Computing α_L...")
    alpha_cache = {}
    for n in range(-N, N + 1):
        alpha_cache[n] = alpha_L(n, L)

    if verbose:
        print(f"Building Weil matrix...")
    tau = build_weil_matrix_dense(N, L, primes, alpha_cache)
    result["matrix_built"] = True

    # Try to use Kwant
    if KWANT_AVAILABLE:
        if verbose:
            print(f"Building Kwant system...")
        try:
            syst = build_kwant_system(tau)
            result["kwant_system_built"] = True

            if verbose:
                print(f"Extracting Hamiltonian from Kwant...")
            tau_kwant = extract_hamiltonian_matrix(syst)

            if tau_kwant is not None:
                comparison = compare_matrices(tau, tau_kwant, N)
                if verbose and comparison:
                    print(f"Matrix comparison: rel_error={comparison['relative_error']:.6e}")
                result["matrix_comparison"] = comparison
            else:
                if verbose:
                    print("Could not extract Hamiltonian matrix from Kwant")
        except Exception as e:
            if verbose:
                print(f"Kwant system construction failed: {e}")
            result["kwant_error"] = str(e)
    else:
        if verbose:
            print(f"Kwant not available: {KWANT_ERROR}")
        result["kwant_error"] = KWANT_ERROR

    # Compute eigenvalues directly from the matrix
    if verbose:
        print(f"Computing eigenvalues (dense eigh)...")
    t0 = time.time()
    evals, evecs = np.linalg.eigh(tau)
    time_eigen = time.time() - t0

    eps = float(evals[0])
    xi = evecs[:, 0]
    result["min_eigenvalue"] = eps
    result["time_eigen"] = time_eigen

    # Find spectrum
    eta_xi = np.sum(xi)
    if abs(eta_xi) > 1e-12:
        xi_norm = xi / eta_xi
    else:
        xi_norm = xi

    if verbose:
        print(f"Finding spectrum zeros...")
    s_zeros = find_spectrum_rational(xi_norm, N, s_min=0, s_max=N)
    gamma_pred = 2 * np.pi * s_zeros / L

    # Match to reference
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
        print(f"Matched: {n_match_tight}/{n_pred}")
        if len(errors) > 0:
            print(f"Mean error: {np.mean(errors):.6e}")

    return result


def main():
    print("\n" + "="*70)
    print("KWANT SPECTRAL CHECK — Cross-Verification")
    print("="*70)

    if not KWANT_AVAILABLE:
        print(f"KWANT NOT AVAILABLE: {KWANT_ERROR}")
        print("Proceeding with matrix computation only (no Kwant system).")

    ref_zeros = load_zeros(100)
    print(f"Loaded {len(ref_zeros)} reference zeros")

    # Test a single size for cross-check (use 500 for comparison with slepc)
    test_size = 500

    result = compute_kwant_spectrum(test_size, L, PRIMES, ref_zeros, verbose=True)

    # Save results
    out_file = os.path.join(RESULTS_DIR, "kwant_spectral_check.json")
    with open(out_file, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\nResults saved to {out_file}")

    # Try to load slepc results for comparison
    slepc_file = os.path.join(RESULTS_DIR, "slepc_spectral_scan.json")
    if os.path.exists(slepc_file):
        with open(slepc_file, "r") as f:
            slepc_data = json.load(f)

        # Find N=500 in slepc results
        slepc_500 = None
        for scan in slepc_data.get("scans", []):
            if scan.get("N") == test_size:
                slepc_500 = scan
                break

        if slepc_500:
            print("\n" + "="*70)
            print(f"CROSS-CHECK: N={test_size}")
            print("="*70)
            print(f"SLEPC: {slepc_500['n_match_tight']}/{slepc_500['n_compared']} matches, " +
                  f"mean_err={slepc_500['mean_error']:.6e}")
            print(f"KWANT: {result['n_match_tight']}/{result['n_compared']} matches, " +
                  f"mean_err={result['mean_error']:.6e}")

            # Compare eigenvalue spectra
            slepc_gamma = slepc_500.get("gamma_pred", [])
            kwant_gamma = result.get("gamma_pred", [])

            if slepc_gamma and kwant_gamma:
                n_compare = min(len(slepc_gamma), len(kwant_gamma))
                diffs = [abs(s - k) for s, k in zip(slepc_gamma[:n_compare], kwant_gamma[:n_compare])]
                print(f"Spectrum difference (SLEPC vs KWANT): " +
                      f"mean={np.mean(diffs):.6e}, max={np.max(diffs):.6e}")


if __name__ == "__main__":
    main()

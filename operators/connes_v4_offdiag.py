"""
Connes v4: Correct OFF-DIAGONAL Archimedean Formula Only

Key finding from v3: the integral formula for the DIAGONAL (2γ_L - 2β_L)
gives wrong values for n≥1 (differs from digamma by ~4). The digamma
diagonal was already correct.

The problem was ONLY the off-diagonal approximation
  (-1)^k · Im[ψ(1/4+is)] / (π|k|)    ← wrong guess, never from paper

Correct off-diagonal from arXiv:2511.22755:
  W_R(V_n, V_m) = [α_L(m) - α_L(n)] / (n - m)

where α_L(n) = (1/π) ∫₀ᴸ sin(2πnx/L) · ρ(x) dx
and   ρ(x)   = e^{x/2} / (2sinh(x))

This file fixes ONLY the off-diagonal while keeping the proven diagonal.
"""

import mpmath
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import json
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR   = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")

mpmath.mp.dps = 80

LAMBDA_SQ = 14
L         = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))
N_TRUNC   = 80
PRIMES    = [2, 3, 5, 7, 11, 13]


def load_zeros(count=100):
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


# ── Archimedean kernel ──────────────────────────────────────────────────────

def rho_scalar(x):
    """ρ(x) = e^{x/2} / (2sinh(x))"""
    if x < 1e-10:
        return 1.0 / (2 * x + 1e-300)
    return np.exp(x / 2) / (np.exp(x) - np.exp(-x))


def alpha_L(n, L):
    """
    α_L(n) = (1/π) ∫₀ᴸ sin(2πnx/L) · ρ(x) dx

    Near x=0: sin(ax) ≈ ax, rho(x) ≈ 1/(2x), so integrand → a/2 = πn/L (finite).
    """
    if n == 0:
        return 0.0
    pi = np.pi
    a  = 2 * pi * n / L

    def integrand(x):
        if x < 1e-12:
            return a / 2.0   # sin(ax)*rho ≈ ax/(2x) = a/2
        return np.sin(a * x) * rho_scalar(x)

    result, _ = quad(integrand, 0, L, limit=300, epsabs=1e-12, epsrel=1e-12,
                     points=[1e-8, 0.01, 0.1, L / 2])
    return result / pi


def build_alpha_cache(L, N):
    """Pre-compute α_L(n) for n = -N to N."""
    cache = {}
    for n in range(-N, N + 1):
        cache[n] = alpha_L(n, L)
    return cache


# ── Diagonal archimedean (digamma formula — proven correct) ─────────────────

def W_arch_diag(n, L):
    """W_R(V_n, V_n) = Re[ψ(1/4 + iπn/L)] + log(2)"""
    if n == 0:
        return float(-mpmath.euler - mpmath.log(2))
    s = float(mpmath.pi) * n / L
    return float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j * s))) + float(mpmath.log(2))


# ── Off-diagonal archimedean (correct formula from paper) ───────────────────

def W_arch_offdiag(n, m, L, alpha_cache):
    """
    W_R(V_n, V_m) = [α_L(m) - α_L(n)] / (n - m)   for n ≠ m
    """
    a_m = alpha_cache[m]
    a_n = alpha_cache[n]
    return (a_m - a_n) / (n - m)


# ── Boundary and Prime terms ─────────────────────────────────────────────────

def W_boundary(n, m, L):
    pi   = float(mpmath.pi)
    L_mp = mpmath.mpf(str(L))
    num  = 32 * L_mp * mpmath.sinh(L_mp / 4)**2 * (L_mp**2 - 16 * pi**2 * m * n)
    den  = (L_mp**2 + 16 * pi**2 * m**2) * (L_mp**2 + 16 * pi**2 * n**2)
    if abs(float(den)) < 1e-50:
        return 0.0
    return float(num / den)


def von_mangoldt(k, primes):
    if k <= 1: return 0.0
    for p in primes:
        if p > k: break
        pk = p
        while pk <= k:
            if pk == k: return float(mpmath.log(p))
            pk *= p
    return 0.0


def W_prime(n, m, L, primes):
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


# ── Build Weil matrix ────────────────────────────────────────────────────────

def build_weil_matrix(N, L, primes, alpha_cache):
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
                arch = W_arch_offdiag(n, m, L, alpha_cache)

            # Boundary
            bnd = W_boundary(n, m, L)

            # Prime
            prm = W_prime(n, m, L, primes)

            val = arch + bnd - prm
            tau[i, j] = val
            tau[j, i] = val

    return tau


# ── ξ̂ zeros ─────────────────────────────────────────────────────────────────

def build_xi_hat_vectorized(z_arr, L, N, alpha_cache):
    """
    ξ̂(z) = (2/√L) sin(zL/2) · Σ_{n=-N}^{N} ξ_n / (z - 2πn/L)
    where ξ_n are eigenvector components.
    """
    pass  # not needed for eigenvalue check


def find_xi_zeros(evec, L, N, z_min=12, z_max=200, n_scan=10000):
    """Find zeros of ξ̂ using the eigenvector of the smallest eigenvalue."""
    pi = np.pi
    ns = np.arange(-N, N + 1)

    def xi_hat(z):
        poles   = 2 * pi * ns / L
        diffs   = z - poles
        safe    = np.where(np.abs(diffs) < 1e-12, 1e-12, diffs)
        weights = evec / safe
        return (2 / np.sqrt(L)) * np.sin(z * L / 2) * np.sum(weights)

    z_scan = np.linspace(z_min, z_max, n_scan)
    vals   = np.array([xi_hat(z) for z in z_scan])

    found_zeros = []
    for k in range(len(vals) - 1):
        if vals[k] * vals[k + 1] < 0:
            try:
                z0 = brentq(xi_hat, z_scan[k], z_scan[k + 1], xtol=1e-10)
                found_zeros.append(z0)
            except Exception:
                pass
    return found_zeros


def main():
    print("CONNES v4: CORRECT OFF-DIAGONAL ARCHIMEDEAN FORMULA")
    print("=" * 60)
    print(f"λ²={LAMBDA_SQ}, L={L:.4f}, N={N_TRUNC}")
    print()

    # Pre-compute α_L values
    print(f"Computing α_L(n) for n = -{N_TRUNC}..{N_TRUNC}...")
    t0 = time.time()
    alpha_cache = build_alpha_cache(L, N_TRUNC)
    print(f"  Done in {time.time()-t0:.1f}s")
    print()

    # Verify diagonal hasn't changed
    print("DIAGONAL VERIFICATION (digamma formula, unchanged):")
    print(f"  {'n':>4}  {'W_arch_diag':>14}")
    for n in range(0, 8):
        print(f"  {n:>4}  {W_arch_diag(n, L):>+14.6f}")
    print()

    # Verify off-diagonal with new formula
    print("OFF-DIAGONAL VERIFICATION (new α_L formula vs old approximation):")
    print(f"  {'(n,m)':>8}  {'Old approx':>12}  {'New α_L':>12}  {'Diff':>10}")
    pi = float(mpmath.pi)
    for n, m in [(1, 0), (2, 0), (3, 1), (5, 2), (-1, 1), (-2, 3)]:
        # Old approximation
        k = n - m
        s_avg = pi * (abs(n) + abs(m)) / (2 * L)
        if s_avg < 0.01: s_avg = 0.01
        psi_im = float(mpmath.im(mpmath.digamma(mpmath.mpf('0.25') + 1j * s_avg)))
        old_val = ((-1)**abs(k)) * psi_im / (pi * abs(k))
        # New formula
        new_val = W_arch_offdiag(n, m, L, alpha_cache)
        print(f"  ({n:>2},{m:>2})  {old_val:>+12.6f}  {new_val:>+12.6f}  {new_val-old_val:>+10.6f}")
    print()

    # Build matrix
    print(f"Building {2*N_TRUNC+1}×{2*N_TRUNC+1} Weil matrix...")
    t0  = time.time()
    tau = build_weil_matrix(N_TRUNC, L, PRIMES, alpha_cache)
    print(f"  Built in {time.time()-t0:.1f}s")

    evals, evecs = np.linalg.eigh(tau)
    eps = float(evals[0])

    print(f"\nε_N = {eps:.8f}")
    print(f"5 smallest: {evals[:5]}")
    print()

    # Compare with previous results
    print("COMPARISON:")
    print(f"  v1  N=30,  diagonal only:      ε_N ≈ +1.144")
    print(f"  v2  N=120, wrong off-diagonal:  ε_N ≈ +0.326")
    print(f"  v4  N={N_TRUNC},  correct off-diag: ε_N = {eps:+.6f}")
    print()

    # Find zeros of ξ̂ using smallest eigenvector
    print("Searching for zeros of ξ̂ in [12, 200]...")
    ref_zeros = load_zeros(100)
    evec_min  = evecs[:, 0]

    found = find_xi_zeros(evec_min, L, N_TRUNC, z_min=12, z_max=200, n_scan=20000)
    print(f"Found {len(found)} zeros")
    print()

    # Match to reference zeros
    tight_tol = 0.5
    loose_tol = 2.0
    n_ref     = min(60, len(ref_zeros))
    n_match_tight = 0
    n_match_loose = 0

    print(f"   #  {'Found':>10}  {'Ref ζ':>10}  {'Error':>10}  Match")
    print("-" * 50)

    for i, ref in enumerate(ref_zeros[:n_ref]):
        if i >= len(found): break
        err = abs(found[i] - ref)
        if err < tight_tol:
            sym = "✓"
            n_match_tight += 1
        elif err < loose_tol:
            sym = "~"
            n_match_loose += 1
        else:
            sym = " "
        if i < 30:
            print(f"  {i+1:>2}  {found[i]:>10.4f}  {ref:>10.4f}  {err:>10.6f}  {sym}")
    n_compared = min(len(found), n_ref)
    print()
    print(f"Matched tight (<{tight_tol}): {n_match_tight}/{n_compared}")
    print(f"Matched loose (<{loose_tol}): {n_match_loose+n_match_tight}/{n_compared}")
    print()

    # ── Plot ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: eigenvalue spectrum
    ax = axes[0]
    ax.plot(evals[:40], 'b^', markersize=6, label=f'Eigenvalues (ε_N={eps:.4f})')
    ax.axhline(0, color='red', linestyle='--', linewidth=1.5, label='ε = 0')
    ax.set_xlabel('Index')
    ax.set_ylabel('Eigenvalue')
    ax.set_title(f'Weil Matrix Spectrum (v4, N={N_TRUNC})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: found zeros vs reference
    ax = axes[1]
    ref_arr   = np.array(ref_zeros[:n_ref])
    found_arr = np.array(found[:n_ref])
    n_plot    = min(len(ref_arr), len(found_arr))
    ax.scatter(found_arr[:n_plot], ref_arr[:n_plot], s=20, c='blue', alpha=0.7,
               label='Predicted vs Reference')
    zmin, zmax = 12, 200
    ax.plot([zmin, zmax], [zmin, zmax], 'r--', linewidth=1, label='Perfect match')
    ax.set_xlabel('Found zero')
    ax.set_ylabel('Reference ζ zero')
    ax.set_title('Zero Matching (v4)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Connes v4: Correct Off-Diagonal | ε_N={eps:.4f} | N={N_TRUNC}', fontsize=11)
    plt.tight_layout()
    out = os.path.join(PLOTS_DIR, "connes_v4_offdiag.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Plot saved: {out}")

    # Save results
    with open(os.path.join(RESULTS_DIR, "connes_v4_offdiag.json"), "w") as f:
        json.dump({
            "N": N_TRUNC, "lambda_sq": LAMBDA_SQ, "L": L,
            "epsilon_N": eps,
            "5_smallest_evals": evals[:5].tolist(),
            "n_found_zeros": len(found),
            "n_match_tight": n_match_tight,
            "n_match_loose": n_match_loose + n_match_tight,
            "n_compared": n_compared,
        }, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()

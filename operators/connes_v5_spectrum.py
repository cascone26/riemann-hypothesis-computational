"""
Connes v5: Correct Spectrum from Lemma 5.4

From arXiv:2511.22755, Lemma 5.4, equation (5.5):

The eigenvalues of the perturbed operator D" (the D_log^{λ,N} in physical notation)
are the zeros s of the rational function:

  P(s) = Σ_{j=-N}^{N}  ξ_j / (j - s)

where ξ_j are the components of the MINIMUM EIGENVECTOR of the Weil matrix τ
(after normalization ⟨η|ξ⟩ = 1, where η = Σ_j V_j).

The physical spectrum (imaginary parts of zeta zeros) is:
  γ_k = (2π/L) · s_k

where s_k are the zeros of P(s).

NOTE: P(s) has poles at integers j = -N,...,N. Between consecutive integers,
there is exactly one zero by continuity (if ξ is "even" = γξ = ξ). So there
are approximately 2N zeros in total. The smallest ~50 should approximate
the first 50 zeta zeros.

Previous approach (wrong ξ̂ formula) gave 1/60 tight matches.
Correct approach should give ~10^{-55} error on first zero, per the paper.
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


# ── Diagonal archimedean (digamma formula) ──────────────────────────────────

def W_arch_diag(n, L):
    """W_R(V_n, V_n) = Re[ψ(1/4 + iπn/L)] + log(2)"""
    if n == 0:
        return float(-mpmath.euler - mpmath.log(2))
    s = float(mpmath.pi) * n / L
    return float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j * s))) + float(mpmath.log(2))


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


# ── Spectrum from Lemma 5.4 ──────────────────────────────────────────────────

def find_spectrum_rational(xi_vec, N, s_min=-80, s_max=80, n_scan=100000):
    """
    Find zeros of P(s) = Σ_{j=-N}^{N} ξ_j / (j - s)

    Poles are at integer j ∈ {-N,...,N}.
    Between consecutive integers, there is (at most) one zero.
    Physical zero = 2π·s/L.
    """
    js = np.arange(-N, N + 1, dtype=float)

    def P(s):
        diffs = js - s
        safe  = np.where(np.abs(diffs) < 1e-12, 1e-12, diffs)
        return float(np.sum(xi_vec / safe))

    # Scan between consecutive poles
    zeros = []
    for j_int in range(int(s_min), int(s_max)):
        # Check in interval (j_int, j_int+1) — but only if both are within our pole range
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


def main():
    print("CONNES v5: CORRECT SPECTRUM VIA LEMMA 5.4")
    print("=" * 60)
    print(f"λ²={LAMBDA_SQ}, L={L:.4f}, N={N_TRUNC}")
    print()

    # Pre-compute α_L values
    print(f"Computing α_L(n) for n = -{N_TRUNC}..{N_TRUNC}...")
    t0 = time.time()
    alpha_cache = {}
    for n in range(-N_TRUNC, N_TRUNC + 1):
        alpha_cache[n] = alpha_L(n, L)
    print(f"  Done in {time.time()-t0:.1f}s")
    print()

    # Build Weil matrix
    print(f"Building {2*N_TRUNC+1}×{2*N_TRUNC+1} Weil matrix...")
    t0  = time.time()
    tau = build_weil_matrix(N_TRUNC, L, PRIMES, alpha_cache)
    print(f"  Built in {time.time()-t0:.1f}s")

    evals, evecs = np.linalg.eigh(tau)
    eps = float(evals[0])
    xi  = evecs[:, 0]   # minimum eigenvector

    print(f"\nε_N = {eps:.8f}")
    print(f"Min eigenvector first 5 components: {xi[:5]}")
    print()

    # Normalize: ⟨η|ξ⟩ = 1 where η = Σ_j V_j → η_j = 1 for all j
    # In discrete basis, ⟨η|ξ⟩ = Σ_j ξ_j
    eta_xi = np.sum(xi)
    if abs(eta_xi) > 1e-12:
        xi_norm = xi / eta_xi
        print(f"Normalized ξ by ⟨η|ξ⟩ = {eta_xi:.6f}")
    else:
        xi_norm = xi
        print(f"WARNING: ⟨η|ξ⟩ = {eta_xi:.2e} (near zero, symmetry issue?)")
    print()

    # Find zeros of P(s) = Σ_j ξ_j/(j-s)
    # Physical zeros are γ = 2π·s/L
    print("Finding spectrum via P(s) = Σ_j ξ_j/(j-s) = 0 ...")
    t0 = time.time()
    js      = np.arange(-N_TRUNC, N_TRUNC + 1, dtype=float)
    s_zeros = find_spectrum_rational(xi_norm, N_TRUNC, s_min=0, s_max=N_TRUNC)
    gamma_pred = 2 * np.pi * s_zeros / L
    print(f"  Found {len(s_zeros)} zeros in s ∈ [0,{N_TRUNC}] in {time.time()-t0:.1f}s")
    print()

    # Load reference zeros
    ref_zeros = load_zeros(100)

    # Match predictions to reference zeros
    tight_tol = 0.5
    loose_tol = 2.0
    n_ref     = min(60, len(ref_zeros))
    n_pred    = min(n_ref, len(gamma_pred))
    n_match_tight = 0
    n_match_loose = 0

    print(f"   #  {'s_zero':>8}  {'Found γ':>10}  {'Ref ζ':>10}  {'Error':>10}  Match")
    print("-" * 60)

    for i in range(min(n_pred, n_ref)):
        ref = ref_zeros[i]
        pred = gamma_pred[i]
        err = abs(pred - ref)
        if err < tight_tol:
            sym = "✓"
            n_match_tight += 1
        elif err < loose_tol:
            sym = "~"
            n_match_loose += 1
        else:
            sym = " "
        if i < 30:
            print(f"  {i+1:>2}  {s_zeros[i]:>8.4f}  {pred:>10.4f}  {ref:>10.4f}  {err:>10.6f}  {sym}")

    print()
    print(f"Matched tight (<{tight_tol}): {n_match_tight}/{n_pred}")
    print(f"Matched loose (<{loose_tol}): {n_match_loose+n_match_tight}/{n_pred}")
    print()

    # Compare with previous approaches
    print("COMPARISON:")
    print(f"  v1  N=30,  diag only:             12/20 (60%)")
    print(f"  v2  N=120, wrong off-diag arch:   28/55 (51%)")
    print(f"  v4  N={N_TRUNC},  correct arch, wrong ξ̂: 1/60 (2%)")
    print(f"  v5  N={N_TRUNC},  correct arch + Lemma5.4: {n_match_tight}/{n_pred} ({100*n_match_tight/max(1,n_pred):.0f}%)")
    print()

    # ── Plot ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: P(s) function
    ax = axes[0]
    s_plot  = np.linspace(0.5, 30.5, 5000)
    P_vals  = []
    for s in s_plot:
        diffs = js - s
        safe  = np.where(np.abs(diffs) < 0.05, np.nan, diffs)
        P_val = np.nansum(xi_norm / safe)
        if np.abs(P_val) > 10: P_val = np.nan
        P_vals.append(P_val)
    P_vals = np.array(P_vals)
    ax.plot(s_plot, P_vals, 'b-', linewidth=0.8, label='P(s)')
    ax.axhline(0, color='red', linewidth=1, linestyle='--')
    for s0 in s_zeros[:20]:
        ax.axvline(s0, color='green', linewidth=0.5, alpha=0.5)
    # Mark reference zero s values
    ref_s = np.array(ref_zeros[:20]) * L / (2 * np.pi)
    for rs in ref_s:
        ax.axvline(rs, color='orange', linewidth=0.5, alpha=0.5, linestyle=':')
    ax.set_xlim(0, 30)
    ax.set_ylim(-5, 5)
    ax.set_xlabel('s (abstract eigenvalue)')
    ax.set_ylabel('P(s) = Σ ξ_j/(j-s)')
    ax.set_title('P(s): green=found, orange=ref ζ zeros (in s units)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: found vs reference
    ax = axes[1]
    n_plot = min(len(gamma_pred), len(ref_zeros), 50)
    ax.scatter(gamma_pred[:n_plot], ref_zeros[:n_plot], s=20, c='blue', alpha=0.7,
               label='Predicted vs Reference')
    zmin, zmax = 12, 200
    ax.plot([zmin, zmax], [zmin, zmax], 'r--', linewidth=1, label='Perfect match')
    ax.set_xlabel('Found zero γ')
    ax.set_ylabel('Reference ζ zero')
    ax.set_title('Zero Matching (v5: Lemma 5.4)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Connes v5: Correct Spectrum | ε_N={eps:.4f} | N={N_TRUNC} | {n_match_tight}/{n_pred} tight', fontsize=11)
    plt.tight_layout()
    out = os.path.join(PLOTS_DIR, "connes_v5_spectrum.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Plot saved: {out}")

    # Save results
    with open(os.path.join(RESULTS_DIR, "connes_v5_spectrum.json"), "w") as f:
        json.dump({
            "N": N_TRUNC, "lambda_sq": LAMBDA_SQ, "L": L,
            "epsilon_N": eps,
            "s_zeros": s_zeros.tolist(),
            "gamma_pred": gamma_pred.tolist(),
            "n_match_tight": n_match_tight,
            "n_match_loose": n_match_loose + n_match_tight,
            "n_compared": n_pred,
        }, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()

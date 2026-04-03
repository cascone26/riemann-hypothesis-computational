"""
Connes v3: Correct Weil Matrix using Exact Archimedean Formula

From arXiv:2511.22755 (Connes-Consani-Moscovici 2025):

OFF-DIAGONAL: W_R(V_n, V_m) = [α_L(m) - α_L(n)] / (n - m)  for n ≠ m

DIAGONAL:     W_R(V_n, V_n) = 2γ_L(n) - 2β_L(n)

where the kernel involves ρ(x) = e^{x/2} / (e^x - e^{-x}) = e^{x/2} / (2sinh(x)):

  α_L(n) = (1/π) ∫₀ᴸ sin(2πnx/L) · ρ(x) dx
  β_L(n) = (1/L) ∫₀ᴸ x · cos(2πnx/L) · ρ(x) dx
  γ_L(n) = ∫₀ᴸ [cos(2πnx/L) - e^{-x/2}] · ρ(x) dx + c(L)

For n=0: α_L(0) = 0 (integrand is 0), W_R(V_0,V_0) = 2γ_L(0) - 2β_L(0)
For n≠0: W_R(V_n,V_n) = 2γ_L(n) - 2β_L(n)
Off-diag: W_R(V_n,V_m) = (α_L(m) - α_L(n)) / (n-m)

The paper reports 10^{-55} error on first zero, 10^{-3} on 50th zero.
Our previous approximation gave ~0.5 error because the off-diagonal was wrong.
"""

import mpmath
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
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
L         = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))  # log(14) ≈ 2.6391
N_TRUNC   = 60   # start smaller to verify correctness

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

def rho(x):
    """ρ(x) = e^{x/2} / (e^x - e^{-x}) = e^{x/2} / (2 sinh x)"""
    if x < 1e-10:
        # ρ(x) → 1/(2x) as x→0
        return 1.0 / (2 * x + 1e-300)
    return np.exp(x/2) / (np.exp(x) - np.exp(-x))


def alpha_L(n, L, epsabs=1e-10):
    """α_L(n) = (1/π) ∫₀ᴸ sin(2πnx/L) ρ(x) dx"""
    if n == 0:
        return 0.0
    pi = np.pi

    def integrand(x):
        if x < 1e-12:
            # sin(2πnx/L)/x → 2πn/L as x→0, so rho*sin ≈ (2πn/L)/(2) = πn/L
            return pi * n / L
        return np.sin(2*pi*n*x/L) * rho(x)

    # Split at potential issues: the integrand has a singularity at x=0
    # but it's integrable (sin(ax)/x → a as x→0 and rho(x) ~ 1/(2x))
    result, err = quad(integrand, 0, L, limit=200, epsabs=epsabs, epsrel=epsabs,
                       points=[1e-8, 0.01, L/2])
    return result / pi


def beta_L(n, L, epsabs=1e-10):
    """β_L(n) = (1/L) ∫₀ᴸ x cos(2πnx/L) ρ(x) dx"""
    pi = np.pi

    def integrand(x):
        if x < 1e-12:
            return 0.5  # x * rho(x) → x/(2x) = 1/2 as x→0, cos → 1
        return x * np.cos(2*pi*n*x/L) * rho(x)

    result, err = quad(integrand, 0, L, limit=200, epsabs=epsabs, epsrel=epsabs,
                       points=[1e-8, 0.01, L/2])
    return result / L


def gamma_L(n, L, c_L=None, epsabs=1e-10):
    """
    γ_L(n) = ∫₀ᴸ [cos(2πnx/L) - e^{-x/2}] ρ(x) dx + c(L)

    The constant c(L) needs to be determined to match W_R(V_0,V_0) = -γ_E - log(2).
    We calibrate it from the n=0 case.
    """
    pi = np.pi

    def integrand(x):
        if x < 1e-12:
            # [cos(2πnx/L) - e^{-x/2}] * rho(x) → [1 - 1] * 1/(2x) → 0/0
            # Need L'Hopital: [cos(ax) - e^{-x/2}] ≈ -a²x²/2 + x/2 for small x
            # times rho(x) ≈ 1/(2x):
            # → (-a²x/2 + 1/2) / 2 → 1/4 as x→0 (for any a)
            a = 2*pi*n/L
            return (1/4.0 - a**2/4.0 * x) if x > 1e-300 else 1/4.0
        cos_part = np.cos(2*pi*n*x/L)
        exp_part = np.exp(-x/2)
        r = rho(x)
        return (cos_part - exp_part) * r

    result, err = quad(integrand, 0, L, limit=200, epsabs=epsabs, epsrel=epsabs,
                       points=[1e-8, 0.01, L/2])

    if c_L is None:
        return result  # return without constant (caller adds it)
    return result + c_L


def compute_arch_correction(L):
    """
    Determine constant c_L by requiring W_R(V_0,V_0) = -γ_E - log(2).
    W_R(V_0,V_0) = 2*gamma_L(0) - 2*beta_L(0) = expected value.
    """
    expected_diag_0 = float(-mpmath.euler - mpmath.log(2))   # ≈ -1.270
    b0 = beta_L(0, L)
    g0_raw = gamma_L(0, L)
    # 2*(g0_raw + c_L) - 2*b0 = expected
    # 2*c_L = expected - 2*g0_raw + 2*b0
    c_L = (expected_diag_0 - 2*g0_raw + 2*b0) / 2.0
    return c_L, g0_raw, b0, expected_diag_0


def W_arch_correct(n, m, L, c_L, alpha_cache=None):
    """
    Correct archimedean Weil term.

    Off-diagonal (n≠m): W_R(V_n, V_m) = [α_L(m) - α_L(n)] / (n - m)
    Diagonal (n=m):     W_R(V_n, V_n) = 2*γ_L(n) - 2*β_L(n)
    """
    if n == m:
        g = gamma_L(n, L) + c_L
        b = beta_L(n, L)
        return 2*g - 2*b
    else:
        if alpha_cache is not None:
            a_m = alpha_cache.get(m)
            a_n = alpha_cache.get(n)
            if a_m is None:
                a_m = alpha_L(m, L)
                alpha_cache[m] = a_m
            if a_n is None:
                a_n = alpha_L(n, L)
                alpha_cache[n] = a_n
        else:
            a_m = alpha_L(m, L)
            a_n = alpha_L(n, L)
        return (a_m - a_n) / (n - m)


# ── Boundary and Prime terms (unchanged from v2) ────────────────────────────

def W_boundary(n, m, L):
    pi = float(mpmath.pi)
    L_mp = mpmath.mpf(str(L))
    num = 32*L_mp * mpmath.sinh(L_mp/4)**2 * (L_mp**2 - 16*pi**2*m*n)
    den = (L_mp**2 + 16*pi**2*m**2) * (L_mp**2 + 16*pi**2*n**2)
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
    pi = float(mpmath.pi)
    k_max = int(float(mpmath.exp(mpmath.mpf(str(L))))) + 1
    total = 0.0
    for k in range(2, k_max + 1):
        lam_k = von_mangoldt(k, primes)
        if lam_k == 0: continue
        y = float(mpmath.log(k))
        weight = lam_k * k**(-0.5)
        if n != m:
            q = (np.sin(2*pi*m*y/L) - np.sin(2*pi*n*y/L)) / (pi*(n-m))
        else:
            q = 2*(1 - abs(y)/L) * np.cos(2*pi*n*y/L)
        total += weight * q
    return total


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    pi = np.pi
    print("CONNES v3: CORRECT ARCHIMEDEAN FORMULA")
    print("=" * 60)
    print(f"λ²={LAMBDA_SQ}, L={L:.4f}, N={N_TRUNC}")
    print()

    # ── 1. Calibrate c_L ──────────────────────────────────────────────
    print("Calibrating archimedean constant c_L...")
    c_L, g0_raw, b0, expected_0 = compute_arch_correction(L)
    print(f"  Expected W_R(V_0,V_0) = {expected_0:.6f}")
    print(f"  gamma_L(0) raw = {g0_raw:.6f}")
    print(f"  beta_L(0)  = {b0:.6f}")
    print(f"  c_L = {c_L:.6f}")
    print(f"  Calibrated W_R(V_0,V_0) = {2*(g0_raw+c_L) - 2*b0:.6f}")
    print()

    # ── 2. Verify diagonal against old formula ────────────────────────
    print("DIAGONAL VERIFICATION (new vs old formula):")
    print(f"{'n':>4}  {'Old (digamma)':>15}  {'New (integral)':>15}  {'Diff':>10}")
    for n_test in range(0, 8):
        if n_test == 0:
            old = float(-mpmath.euler - mpmath.log(2))
        else:
            s = pi * n_test / L
            old = float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j*s))) + float(mpmath.log(2))
        new = W_arch_correct(n_test, n_test, L, c_L)
        print(f"{n_test:>4}  {old:>15.6f}  {new:>15.6f}  {new-old:>10.6f}")
    print()

    # ── 3. Build Weil matrix ─────────────────────────────────────────
    print(f"Building {2*N_TRUNC+1}×{2*N_TRUNC+1} Weil matrix...")
    dim = 2*N_TRUNC + 1
    tau = np.zeros((dim, dim))
    alpha_cache = {}

    t0 = time.time()
    for i in range(dim):
        n = i - N_TRUNC
        for j in range(i, dim):
            m = j - N_TRUNC
            arch = W_arch_correct(n, m, L, c_L, alpha_cache)
            bnd  = W_boundary(n, m, L)
            prm  = W_prime(n, m, L, PRIMES)
            val  = arch + bnd - prm
            tau[i, j] = val
            tau[j, i] = val
    print(f"  Built in {time.time()-t0:.1f}s")
    print()

    # ── 4. Eigendecomposition ────────────────────────────────────────
    evals, evecs = np.linalg.eigh(tau)
    eps = float(evals[0])
    xi  = evecs[:, 0].copy()
    delta = xi.sum() / np.sqrt(L)
    if abs(delta) > 1e-12:
        xi /= delta

    print(f"ε_N = {eps:.8f}")
    print(f"5 smallest: {evals[:5]}")
    print()

    # Compare to v2 (which used approx off-diagonal):
    # v2 had ε_N ≈ 0.326 for N=120. Correct formula should give much smaller ε_N.

    # ── 5. Find zeros of ξ̂ ────────────────────────────────────────
    poles = np.array([2*pi*j/L for j in range(-N_TRUNC, N_TRUNC+1)])

    def xi_hat(z):
        prefix = (2.0/np.sqrt(L)) * np.sin(z*L/2.0)
        diffs  = z - poles
        mask   = np.abs(diffs) > 1e-12
        return prefix * np.sum(xi[mask] / diffs[mask])

    print("Searching for zeros of ξ̂ in [12, 200]...")
    ref_zeros = load_zeros(100)
    z_grid = np.linspace(12.0, 200.0, 80000)
    vals   = np.array([xi_hat(z) for z in z_grid])

    found = []
    for i in range(len(z_grid)-1):
        if vals[i]*vals[i+1] < 0 and not (np.isnan(vals[i]) or np.isnan(vals[i+1])):
            try:
                root = brentq(xi_hat, z_grid[i], z_grid[i+1], xtol=1e-9)
                found.append(root)
            except Exception:
                pass

    print(f"Found {len(found)} zeros")
    print()

    # Match against reference
    matched_tight = 0
    matched_loose = 0
    print(f"{'#':>4}  {'Found':>10}  {'Ref ζ':>10}  {'Error':>10}  {'Match':>6}")
    print("-" * 50)
    for i, fz in enumerate(found[:60]):
        diffs = np.abs(ref_zeros - fz)
        idx   = np.argmin(diffs)
        err   = diffs[idx]
        tag   = "✓" if err < 0.5 else ("~" if err < 2.0 else " ")
        if err < 0.5: matched_tight += 1
        if err < 2.0: matched_loose += 1
        if i < 30:
            print(f"{i+1:>4}  {fz:>10.4f}  {ref_zeros[idx]:>10.4f}  {err:>10.6f}  {tag:>6}")

    n_shown = min(len(found), 60)
    print()
    print(f"Matched tight (<0.5): {matched_tight}/{n_shown}")
    print(f"Matched loose (<2.0): {matched_loose}/{n_shown}")
    print()
    print("Previous results:")
    print("  v1  N=30,  diag only:     12/20 (60%)")
    print("  v2  N=120, wrong off-diag: 28/55 (51%)")
    print(f"  v3  N={N_TRUNC}, correct:     {matched_tight}/{n_shown} ({100*matched_tight/max(n_shown,1):.0f}%)")

    # Save results
    results = {
        "N": N_TRUNC, "lambda_sq": LAMBDA_SQ, "L": L, "c_L": c_L,
        "epsilon_N": eps, "eigenvalues_5": evals[:5].tolist(),
        "n_found": len(found), "matched_tight": matched_tight,
        "matched_loose": matched_loose, "n_shown": n_shown,
    }
    with open(os.path.join(RESULTS_DIR, "connes_v3.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()

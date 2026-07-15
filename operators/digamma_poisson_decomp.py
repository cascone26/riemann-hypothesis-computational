#!/usr/bin/env python3
"""
digamma_poisson_decomp.py — Explore the Poisson decomposition of Omega(r).

The key object in A(g) = (1/2pi) int Omega(r)|g_hat(1/2+ir)|^2 dr is:

  Omega(r) = Re psi(1/4 + ir/2) - log(pi)

By the partial fraction / Mittag-Leffler formula for psi:

  psi(s) = -gamma - 1/s + sum_{n=1}^inf (1/n - 1/(n+s))

So:
  Re psi(1/4 + ir/2) = Re psi(1/4 + ir/2)

Using the digamma reflection formula and product formula:
  Omega(r) = sum_{n=0}^inf [ (n+1/4) / ((n+1/4)^2 + r^2/4) ] - log(pi) - ???

Actually from the Hadamard product formula for xi(s):
  log xi(s) = (rank terms) + sum_rho log(1 - s/rho)

And from the relation between Omega and the explicit formula:
  Omega(r) = -gamma - log(pi) + sum_{n=0}^inf [ 1/(n+1/4+ir/2) + 1/(n+1/4-ir/2) - 2/n ] / 2

Let me verify this numerically. The key formula is:
  psi(s) = -gamma + sum_{n=0}^inf (1/(n+1) - 1/(n+s))

So:
  psi(1/4 + ir/2) = -gamma + sum_{n=0}^inf (1/(n+1) - 1/(n+1/4+ir/2))

  Re psi(1/4 + ir/2) = -gamma + sum_{n=0}^inf (1/(n+1) - Re[1/(n+1/4+ir/2)])
  = -gamma + sum_{n=0}^inf (1/(n+1) - (n+1/4)/((n+1/4)^2 + r^2/4))

  Omega(r) = Re psi(1/4+ir/2) - log pi
  = -gamma - log pi + sum_{n=0}^inf (1/(n+1) - (n+1/4)/((n+1/4)^2+r^2/4))

The term (n+1/4)/((n+1/4)^2+r^2/4) is a POSITIVE Lorentzian (Poisson kernel type).

So:
  Omega(r) = C_0 + sum_{n=0}^inf [ 1/(n+1) - P_n(r) ]

where C_0 = -gamma - log pi and P_n(r) = (n+1/4)/((n+1/4)^2+r^2/4) > 0.

At r=0: P_n(0) = 1/(n+1/4), so:
  Omega(0) = C_0 + sum_{n=0}^inf (1/(n+1) - 1/(n+1/4))
  = C_0 + sum_{n=0}^inf [ 1/(n+1) - 1/(n+1/4) ]  (this series converges)
  = C_0 + [psi(1) - psi(1/4)] = (-gamma-logpi) + (-gamma) - (-gamma - pi/2 - 3*log2)
  = -logpi + pi/2 + 3log2   [verify numerically]

The integral:
  int P_n(r) |g_hat(1/2+ir)|^2 dr = int (n+1/4)/((n+1/4)^2+r^2/4) |g_hat|^2 dr

This is a WEIGHTED AVERAGE of |g_hat|^2 with a Lorentzian weight. This integral is POSITIVE.

So A(g) decomposes as:
  A(g) = (1/2pi) * C_0 * ||g_hat||^2 + (1/2pi) * sum_n [ (1/(n+1)) * ||g_hat||^2 - int P_n|g_hat|^2 ]

The first term C_0 * ||g_hat||^2 is NEGATIVE (C_0 = -gamma-logpi < 0).
The Poisson terms can be positive or negative depending on how g_hat varies.

KEY QUESTION: Can we bound each Lorentzian integral from below?
  int P_n(r) |g_hat(1/2+ir)|^2 dr ≤ (1/(n+1/4)) * int_{-infty}^infty |g_hat|^2 dr (at r=0 peak)
  But also ≥ something related to the mass of g_hat.

This script computes:
1. The Poisson decomposition of Omega(r) numerically
2. The contribution of each Lorentzian term to A(g) for a specific g
3. Whether there's a pattern of cancellation between terms
"""

import json, os
import numpy as np
from scipy.special import digamma as _digamma
from scipy.integrate import quad
import mpmath as mp

mp.mp.dps = 50

LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
GAMMA = float(mp.euler)
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/digamma_poisson_decomp.json")


# === Verify the Poisson decomposition formula for Omega ===

def omega_exact(r):
    """Omega(r) = Re psi(1/4+ir/2) - log pi"""
    if hasattr(r, '__len__'):
        return np.real(_digamma(0.25 + 0.5j * np.array(r))) - LOGPI
    return float(np.real(_digamma(0.25 + 0.5j * r))) - LOGPI


def omega_poisson_expansion(r, N_terms=200):
    """
    Omega(r) = -gamma - logpi + sum_{n=0}^{N_terms} (1/(n+1) - (n+1/4)/((n+1/4)^2+r^2/4))
    + tail correction
    """
    r2 = r**2 / 4
    total = -GAMMA - LOGPI
    for n in range(N_terms):
        a = n + 0.25
        total += 1.0/(n+1) - a/(a**2 + r2)
    return total


def poisson_term(r, n):
    """P_n(r) = (n+1/4) / ((n+1/4)^2 + r^2/4)  (Lorentzian)"""
    a = n + 0.25
    if hasattr(r, '__len__'):
        r = np.array(r)
        return a / (a**2 + r**2/4)
    return a / (a**2 + r**2/4)


def g_hat_constant(r, norm_L=True):
    """
    g = 1/sqrt(L) * rect on [0,L] (roughly the n=0 mode)
    g_hat(1/2+ir) ~ sqrt(L) * sinc(rL/2) = sqrt(L) * sin(rL/2)/(rL/2)

    Actually in the Connes basis, n=0 mode ehat_0 is:
    ehat_0(r) = sqrt(L) when r=0, (1/sqrt(L)) * 2sin(rL/2)/r elsewhere

    For simplicity use ehat_0 directly as the Fourier transform.
    """
    eps = 1e-10
    if hasattr(r, '__len__'):
        r = np.array(r)
        out = np.empty_like(r, dtype=float)
        near = np.abs(r) < eps
        out[near] = np.sqrt(L)
        d = r[~near]
        out[~near] = (2.0/np.sqrt(L)) * np.sin(d * L/2) / d
        return out
    if abs(r) < eps:
        return np.sqrt(L)
    return (2.0/np.sqrt(L)) * np.sin(r * L/2) / r


# === Main analysis ===

def analyze_decomposition():
    results = {}

    # 1. Verify Poisson decomposition at several r values
    print("=== 1. POISSON DECOMPOSITION VERIFICATION ===", flush=True)
    r_probe = [0.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    decomp_check = []
    for r in r_probe:
        exact = omega_exact(r)
        approx = omega_poisson_expansion(r, N_terms=500)
        decomp_check.append({"r": r, "exact": exact, "poisson_approx": approx, "diff": exact - approx})
        print(f"  r={r:6.1f}: exact={exact:.8f}  poisson_approx={approx:.8f}  diff={exact-approx:.2e}", flush=True)
    results["decomp_verification"] = decomp_check

    # 2. For g = ehat_0 (n=0 mode), compute contribution of each Lorentzian to A(g)
    print("\n=== 2. LORENTZIAN CONTRIBUTION TO A(ehat_0) ===", flush=True)

    # A(g) = (1/2pi) int Omega(r) |g_hat(r)|^2 dr
    # Decompose Omega using Poisson series:
    # A(g) = (1/2pi) * [ (-gamma-logpi) int |g|^2 + sum_n int (1/(n+1)) |g|^2 - sum_n int P_n |g|^2 ]

    def integrand_full(r):
        return omega_exact(r) * g_hat_constant(r)**2 / (2*np.pi)

    # Compute A(ehat_0) directly
    A_direct, err = quad(integrand_full, 0, 200, limit=300, epsabs=1e-12)
    A_direct *= 2  # symmetrize (integral is even)
    print(f"  A(ehat_0) direct = {A_direct:.8f} (err={err:.2e})", flush=True)

    # Compute ||ehat_0||^2 (the L2 norm)
    def norm_sq(r):
        return g_hat_constant(r)**2
    norm_sq_val, _ = quad(norm_sq, 0, 200, limit=300)
    norm_sq_val *= 2
    print(f"  ||ehat_0||^2 = {norm_sq_val:.6f} (should be L={L:.4f})", flush=True)

    # Constant term contribution: C_0 * ||g||^2 / (2pi)
    C0 = -GAMMA - LOGPI
    C0_contribution = C0 * norm_sq_val / (2*np.pi)
    print(f"  C_0 = -gamma - log(pi) = {C0:.8f}", flush=True)
    print(f"  C_0 contribution to A = {C0_contribution:.8f}", flush=True)

    # Lorentzian contributions: for each n, compute (1/(n+1) - int P_n |g|^2) * ||g||^2...
    # Actually: contribution_n = (1/2pi) int (1/(n+1) - P_n(r)) |g_hat(r)|^2 dr
    print(f"\n  Term-by-term Lorentzian contributions (first 20 terms):", flush=True)
    lorentz_terms = []
    cumulative = C0_contribution
    for n in range(20):
        # int P_n(r) |g_hat(r)|^2 / (2pi) dr
        def P_integrand(r, n=n):
            return poisson_term(r, n) * g_hat_constant(r)**2 / (2*np.pi)
        P_int, _ = quad(P_integrand, 0, 200, limit=300)
        P_int *= 2

        # term contribution = (1/(n+1)) * ||g||^2 / (2pi) - P_int
        const_n = norm_sq_val / ((n+1) * 2 * np.pi)
        term = const_n - P_int
        cumulative += term
        lorentz_terms.append({"n": n, "const_part": const_n, "P_int": P_int, "term": term, "cumulative": cumulative})
        print(f"    n={n:2d}: const={const_n:.6f}  P_int={P_int:.6f}  term={term:+.6f}  cumul={cumulative:.6f}", flush=True)

    results["A_ehat0_direct"] = A_direct
    results["norm_sq_ehat0"] = norm_sq_val
    results["C0"] = C0
    results["C0_contribution"] = C0_contribution
    results["lorentzian_terms"] = lorentz_terms

    # 3. Compute B(ehat_0) and P(ehat_0) for comparison
    # B(ehat_0) = W_{0,2}[0,0] (the (0,0) element of the boundary term)
    sh = np.sinh(L/4)**2
    B_00 = 32 * L * sh * (L**2) / L**4  # = 32*sinh^2(L/4)/L
    print(f"\n=== 3. B(ehat_0) and P(ehat_0) ===", flush=True)
    print(f"  B(ehat_0) = B[0,0] = {B_00:.8f}", flush=True)

    log_k = np.array([np.log(k) for k, _ in [(k, lk) for k in range(2, LAMBDA_SQ+1) for lk in [0] if (lambda k: any(k == p**m for p in [2,3,5,7,11,13] for m in range(1,20) if p**m <= k and k == p**m))(k)]])
    # simpler: rebuild VM
    VM = []
    for k in range(2, LAMBDA_SQ+1):
        n, p, f, ok = k, 0, 2, True
        while f*f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            VM.append((k, float(mp.log(p))))

    P_00 = sum(lk * k**(-0.5) * 2 * (1 - np.log(k)/L) for k, lk in VM)
    print(f"  P(ehat_0) = P[0,0] = {P_00:.8f}", flush=True)
    print(f"  A+B-P for ehat_0 = {A_direct + B_00 - P_00:.8e}", flush=True)
    print(f"  (This is tau[0,0] from the finite basis formula at N->inf)", flush=True)

    results["B_00"] = B_00
    results["P_00"] = P_00
    results["tau_00_analytic"] = A_direct + B_00 - P_00

    # 4. Check: is Omega(r) "mostly positive" for the right test function?
    print("\n=== 4. OMEGA SIGN ANALYSIS ===", flush=True)
    r_fine = np.linspace(0, 200, 50000)
    Om = omega_exact(r_fine)
    zero_crossing = float(r_fine[np.argmin(np.abs(Om))])
    print(f"  Omega zero crossing at r ≈ {zero_crossing:.4f}", flush=True)
    print(f"  Omega(0) = {omega_exact(0):.6f}", flush=True)
    print(f"  Omega(100) = {omega_exact(100):.6f}", flush=True)
    print(f"  Omega at large r ≈ log(r/2) - gamma (grows to +inf)", flush=True)

    with open(OUTFILE, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults: {OUTFILE}", flush=True)
    return results


def main():
    print("DIGAMMA POISSON DECOMPOSITION ANALYSIS", flush=True)
    print(f"lambda^2={LAMBDA_SQ}, L={L:.4f}", flush=True)
    analyze_decomposition()


if __name__ == "__main__":
    main()

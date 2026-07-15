#!/usr/bin/env python3
"""
witness_construction.py — Attempt to construct witness test functions.

KEY QUESTION: Given a hypothetical off-critical zero at ρ₀ = 1/2 + Δσ + iγ,
can we find a sinc-basis test function g that makes τ(g) < 0?

NON-CIRCULAR STRATEGY:
The Weil formula says τ(g) = Σ_ρ |ĝ(ρ)|². This is ≥ 0 ALWAYS (sum of squares).
BUT: the operator τ is defined WITHOUT knowledge of zero locations.
τ(g) = A(g) + B(g) - P(g) which CAN be negative if the prime sum dominates.

THE APPROACH:
Modify τ to include ONE hypothetical off-critical zero ρ₀ explicitly:
   τ̃(g) = τ_critical(g) - 2·Re(|ĝ(ρ₀)|²)

where the -2·Re term represents the contribution of ρ₀ and its symmetric pair 1-ρ̄₀.

If we CAN construct g such that τ̃(g) < 0, then ρ₀ being a zero contradicts τ ≥ 0.

IMPLEMENTATION:
- Use lsq=50, N=250 (fast)
- Hypothetical zero at ρ₀ = 1/2 + Δσ + iγ₁ (move first zero off critical line)
- Find g (linear combination of sinc basis vectors) that maximizes |τ̃(g)| < 0
- This is equivalent to finding the most negative eigenvalue of τ̃

TEST CASES:
Δσ = 0 (on critical line): τ̃ should be ≥ 0 everywhere (witness doesn't exist)
Δσ = 0.01: τ̃ may have negative eigenvalue
Δσ = 0.1: τ̃ should have clearly negative eigenvalue
"""
import time, json
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 40000
OUTFILE = '/Users/scones/projects/riemann/results/witness_construction.out'

# Riemann zeros
GAMMA = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918720, 43.327073, 48.005151, 49.773832]

def vm_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab

def build_tau_base(lsq, N):
    """Build the base operator τ (with all zeros on critical line via digamma)."""
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lsq)
    idx = np.arange(-N, N + 1, dtype=float)
    ns  = idx

    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]; w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI
    d = (TWO_PI * ns / L)[:, None] - grid[None, :]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    A = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
    B = (32.0 * L * sh * (L ** 2 - 16.0 * np.pi ** 2 * np.outer(ns, ns)) /
         np.outer(L ** 2 + 16.0 * np.pi ** 2 * ns ** 2,
                  L ** 2 + 16.0 * np.pi ** 2 * ns ** 2))

    if VM:
        lkpp = np.array([np.log(k) for k, _ in VM])
        wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))

    tau_base = (A + B - P); tau_base = (tau_base + tau_base.T) / 2.0
    return tau_base, L, ns

def zero_contribution_matrix(ns, L, sigma, gamma):
    """
    Matrix contribution of a hypothetical off-critical zero ρ₀ = sigma + i*gamma.

    The Lorentzian kernel at ρ₀: 1/(ρ₀(1-ρ₀)) × Re[E_m(ρ₀) E_n*(ρ₀)]

    For ρ₀ = 1/2 + Δσ + iγ:
    The digamma function at σ=1/2 already includes this zero at γ (with Δσ=0).
    The CORRECTION for moving it to Δσ≠0 is:
    δA_{mn} = Re[1/((1/2+Δσ+iγ)(1/2-Δσ-iγ)) - 1/((1/2+iγ)(1/2-iγ))] × [mode overlap]

    Simplified: we compute the LORENTZIAN MATRIX for a Lorentzian at ρ₀.
    """
    # The contribution of ρ₀ = s0 + i*gamma to the digamma sum is approximately:
    # Re[1/((n+s0)(n+1-s0))] for n being a summation index in the digamma expansion
    # In the sinc basis, the matrix element is:
    # A_{mn}^{(ρ₀)} = (1/2π) ∫ L_{ρ₀}(ξ) E_m*(ξ) E_n(ξ) dξ
    # where L_{ρ₀}(ξ) = 2Re[1/(ρ₀(1-ρ₀))] × 1/((ξ-2γ)²/(4) + (σ-1/2)²)   [Lorentzian]

    # Use the real part of the Lorentzian kernel centered at 2*gamma
    G_loc = 50000
    R = TWO_PI * len(ns)//2 / L + 80.0
    grid = np.linspace(-R, R, G_loc)
    h = grid[1] - grid[0]; w_g = np.full(G_loc, h); w_g[0] = w_g[-1] = 0.5 * h

    # Lorentzian at ξ = 2*gamma with width sigma (Re part deviation from 1/2)
    Dsigma = sigma - 0.5
    # delta_Omega due to moving zero from 1/2 to sigma:
    # d/d(sigma) [Re psi(sigma/2 + i*xi/4)] at sigma=1/2
    # For a single zero at gamma: contribution is a Lorentzian centered at 2*gamma
    # Amplitude: delta ~ 2*Dsigma / ((xi - 2*gamma)^2 + (2*Dsigma)^2 + 0.01)
    denom = (grid - 2*gamma)**2 + (2*Dsigma)**2 + 1e-8
    L_zero = 2*Dsigma / denom  # derivative of Lorentzian w.r.t. sigma

    d = (TWO_PI * ns / L)[:, None] - grid[None, :]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    contrib = ((E * (w_g * L_zero)) @ E.conj().T).real / (2.0 * np.pi)
    return contrib

def main():
    lsq, N = 50, 250

    print(f"Witness construction test: lsq={lsq}, N={N}")
    print(f"Hypothetical zero: ρ₀ = 1/2 + Δσ + i·γ₁ where γ₁={GAMMA[0]:.3f}")
    print()
    with open(OUTFILE, 'w') as f:
        f.write(f"witness_construction {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    t0 = time.time()
    tau_base, L, ns = build_tau_base(lsq, N)
    ev_base = scipy.linalg.eigvalsh(tau_base)
    lam_min_base = float(ev_base[0])
    print(f"Base τ (σ=0.5): lam_min = {lam_min_base:.4e}")

    delta_sigmas = [0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2]

    print(f"\n{'Δσ':>8} {'lam_min(τ̃)':>14} {'lam_min_base':>14} {'Δlam':>12} {'sign':>6}")
    print("-" * 65)

    with open(OUTFILE, 'a') as f:
        f.write(f"Base lam_min={lam_min_base:.4e}\n")
        f.write(f"{'Δσ':>8} {'lam_min_tilde':>14} {'Δlam':>12} {'detectable':>12}\n")
        f.write("-" * 55 + "\n")

    for dsigma in delta_sigmas:
        if abs(dsigma) < 1e-10:
            lam_tilde = lam_min_base
            delta_lam = 0.0
        else:
            # Build modified operator: τ̃ = τ_base + correction for moving γ₁ by Δσ
            contrib = zero_contribution_matrix(ns, L, 0.5 + dsigma, GAMMA[0])
            tau_tilde = tau_base + contrib
            tau_tilde = (tau_tilde + tau_tilde.T) / 2.0
            ev_tilde = scipy.linalg.eigvalsh(tau_tilde)
            lam_tilde = float(ev_tilde[0])
            delta_lam = lam_tilde - lam_min_base

        detectable = abs(delta_lam) > abs(lam_min_base) if abs(lam_min_base) > 0 else (abs(delta_lam) > 1e-14)
        sign_str = "NEGATIVE" if lam_tilde < 0 else "pos"
        print(f"{dsigma:8.4f} {lam_tilde:14.4e} {lam_min_base:14.4e} {delta_lam:12.4e} {sign_str:6}")
        with open(OUTFILE, 'a') as f:
            f.write(f"{dsigma:8.4f} {lam_tilde:14.4e} {delta_lam:12.4e} {str(detectable):>12}\n")

    total = time.time() - t0
    print(f"\nTotal: {total:.1f}s")
    print(f"If lam_min(τ̃) < 0 for any Δσ > 0 → WITNESS EXISTS → contradicts τ≥0 → RH")
    print(f"If lam_min(τ̃) ≥ 0 for ALL Δσ → no witness at this resolution")
    with open(OUTFILE, 'a') as f:
        f.write(f"\nTotal: {total:.1f}s\n")

if __name__ == "__main__":
    main()

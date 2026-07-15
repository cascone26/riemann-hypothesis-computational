#!/usr/bin/env python3
"""
analytic_slope_verify.py — Verify dlam_min/dσ analytically.

From perturbation theory:
  dlam_min/dσ = ⟨ψ₀, dτ/dσ ψ₀⟩ = -⟨ψ₀, dA/dσ ψ₀⟩

where dA/dσ is the derivative of the A matrix w.r.t. the σ-position of zero γ₁.

dA_{mn}/dσ|_{σ=1/2} = -∫ G(ξ,γ₁) × E_m(ξ) E_n^*(ξ) dξ / (2π) + [conjugate from rho_bar]

where G(ξ,γ₁) = d/dσ Re[1/(1/4+iξ/2 - σ - iγ₁)]|_{σ=1/2}
             = -Re[d/dσ 1/(σ-1/4-i(ξ/2-γ₁))]|_{σ=1/2}
             = Re[1/(1/4-σ + i(ξ/2-γ₁))²]|_{σ=1/2}
             = Re[1/(i(ξ/2-γ₁) - (-1/4))²]  -- wait let me redo

At σ_0=1/2: denominator s = σ₀ - 1/4 + i(ξ/2-γ₁) = 1/4 + i(ξ/2-γ₁)
            1/s = (1/4 - i(ξ/2-γ₁)) / ((1/4)² + (ξ/2-γ₁)²)
d/dσ (1/s) = -1/s²
So d/dσ Re[1/s]|_{σ=1/2} = Re[-1/s²] = -Re[1/s²]

s² = (1/4)² - (ξ/2-γ₁)² + 2i(1/4)(ξ/2-γ₁)
1/s² = ((1/4)² - (ξ/2-γ₁)² - 2i(1/4)(ξ/2-γ₁)) / ((1/4)²+(ξ/2-γ₁)²)²

Re[1/s²] = ((1/4)²-(ξ/2-γ₁)²) / ((1/4)²+(ξ/2-γ₁)²)²

So:
G(ξ,γ₁) = -Re[1/s²] = -((1/4)²-(ξ/2-γ₁)²) / ((1/4)²+(ξ/2-γ₁)²)²
         = ((ξ/2-γ₁)²-(1/4)²) / ((1/4)²+(ξ/2-γ₁)²)²

This is POSITIVE for |ξ/2-γ₁| > 1/4 and NEGATIVE for |ξ/2-γ₁| < 1/4.
Near the zero frequency ξ≈2γ₁: G is NEGATIVE (the integrand is negative near the zero).

So dA_{mn}/dσ = -∫ [G_fwd + G_bwd] × E_m E_n^* dξ / (2π)
where G_fwd from ρ₀ = 1/2+iγ₁ and G_bwd from ρ₀_bar = 1/2-iγ₁.

The KEY question: is ⟨ψ₀, dA/dσ ψ₀⟩ > 0?
If yes → d/dσ lam_min = -⟨ψ₀, dA/dσ ψ₀⟩ < 0 → slope NEGATIVE (confirmed by numerics).
"""
import os, json, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
THRESH = 0.02
GAMMA1 = 14.134725
OUTFILE = os.path.expanduser("~/riemann/results/analytic_slope_verify.json")


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


def build_and_analyze(lsq, N):
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
    arch = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
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

    AB = arch + B
    tau = (AB - P); tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    lam0 = float(ev[0])
    # Ground state (minimum eigenvalue eigenvector)
    psi0 = V[:, 0]  # 1D vector, already normalized

    # Compute dA/dσ analytically at σ₀=1/2 for γ=GAMMA1
    # G(ξ) = ((ξ/2-γ₁)² - (1/4)²) / ((1/4)² + (ξ/2-γ₁)²)²  [for ρ₀=1/2+iγ₁]
    # Plus conjugate from ρ₀_bar = 1/2-iγ₁ (replace γ₁ → -γ₁):
    # G_bar(ξ) = ((ξ/2+γ₁)² - (1/4)²) / ((1/4)² + (ξ/2+γ₁)²)²
    u_fwd = grid / 2 - GAMMA1
    u_bwd = grid / 2 + GAMMA1
    quarter_sq = (0.25) ** 2
    G_fwd = (u_fwd ** 2 - quarter_sq) / (quarter_sq + u_fwd ** 2) ** 2
    G_bwd = (u_bwd ** 2 - quarter_sq) / (quarter_sq + u_bwd ** 2) ** 2
    dOm_dsigma = G_fwd + G_bwd  # total dΩ/dσ at σ=1/2

    # dA/dσ matrix: ∫ dΩ/dσ × E_m E_n^* dξ / (2π)
    dA_dsigma = ((E * (w * dOm_dsigma)) @ E.conj().T).real / (2.0 * np.pi)

    # Analytic slope = -⟨ψ₀, dA/dσ ψ₀⟩
    analytic_slope = -float(psi0 @ dA_dsigma @ psi0)

    # Also project onto near-null subspace
    null_mask = np.abs(ev) < THRESH
    dim = null_mask.sum()
    Vn = V[:, null_mask]  # near-null basis
    # Project dA onto near-null: ⟨ψ_null, dA ψ_null⟩ expectation in null subspace
    dA_null_proj = float(np.trace(Vn.T @ dA_dsigma @ Vn)) / max(1, dim)
    analytic_slope_null = -dA_null_proj

    # Numerical slope (from finite difference)
    ds = 0.001
    s0_g = 0.25 + 0.5j * grid
    def contrib(sigma):
        rho = sigma + 1j * GAMMA1
        rho_bar = sigma - 1j * GAMMA1
        return (np.real(1.0 / (s0_g - rho) + 1.0 / rho) +
                np.real(1.0 / (s0_g - rho_bar) + 1.0 / rho_bar))
    dOm_fd = contrib(0.5 + ds) - contrib(0.5)
    dA_fd = ((E * (w * dOm_fd)) @ E.conj().T).real / (2.0 * np.pi)
    tau2 = (AB + dA_fd - P); tau2 = (tau2 + tau2.T) / 2.0
    lam_g = float(scipy.linalg.eigvalsh(tau2, subset_by_index=[0,0])[0])
    numeric_slope = (lam_g - lam0) / ds

    # Agreement check
    agree_ratio = analytic_slope / numeric_slope if numeric_slope != 0 else float('nan')

    return {
        'lsq': lsq, 'N': N, 'L': round(L, 5),
        'lam0': lam0, 'dim': int(dim),
        'analytic_slope': float(analytic_slope),
        'analytic_slope_null': float(analytic_slope_null),
        'numeric_slope': float(numeric_slope),
        'agree_ratio': float(agree_ratio),
    }


def main():
    lsq_list = [5, 8, 12, 20, 30, 50, 80, 100]
    results = []
    t0 = time.time()
    print(f"{'lsq':>5}  {'analytic':>12}  {'numeric':>12}  {'ratio':>8}  {'match?':>8}")
    for lsq in lsq_list:
        N = max(20, 5 * lsq)
        t1 = time.time()
        r = build_and_analyze(lsq, N)
        dt = time.time() - t1
        print(f"lsq={r['lsq']:4d}: analytic={r['analytic_slope']:+.5e}  numeric={r['numeric_slope']:+.5e}  "
              f"ratio={r['agree_ratio']:.4f}  ({dt:.1f}s)", flush=True)
        r['dt'] = round(dt, 1)
        results.append(r)

    out = {'gamma1': GAMMA1, 'delta_sigma': 0.001, 'results': results,
           'total_time': round(time.time() - t0, 1)}
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Investigate whether τ = W02 - WR - Wp belongs to a known universality class
from random matrix theory, particularly Toeplitz-plus-low-rank matrices.

Key structure of τ:
1. W02 and WR are dense, smooth kernels (Toeplitz-like structure from Fourier basis)
2. Wp is a sum of M = O(c/ln(c)) rank-1 or rank-2 perturbations (from prime contributions)
3. Negative eigenvalues at t=0 (pure background) decrease as primes are added (t→1)

Classical results on extreme eigenvalues of structured matrices:

A. **Toeplitz matrices:** if T is a Toeplitz matrix with symbol f(θ), then
   λ_min(T) ~ (2π)^(-1) min_{θ∈[0,2π]} f(θ)

B. **Toeplitz + finite-rank perturbation:** If T is Toeplitz and P is rank-r, then
   the smallest eigenvalue of T+P depends on the angle between the rank-r perturbation
   and the near-null space of T.

C. **Perturbation theory for nearly-singular perturbations:** If λ is small and
   v is the corresponding eigenvector, then by 1st-order perturbation:
   δλ ≈ ⟨v, ΔT v⟩ where ΔT is the perturbation.

Our case:
- W02-WR at t=0 has negative eigenvalues (spectrum in negative half-line)
- Wp as t increases adds positive-definite structure (primes contribute positively)
- At t=1, the smallest eigenvalue is tiny but positive

The mechanism: the "near-null space" at t=0 (small |λ| negative eigenvalues) gets
perturbedby the primes, pushing them to near-zero positive values.

Key question: **How does the minimum eigenvalue scale with the rank (M) and
magnitude of Wp vs the background negativity?**

Hypothesis (to test):
  λ_min ~ (background-negativity) / (prime-cancellation-effectiveness)^α

where α might be related to the spectral norm σ(Wp) and the gap structure.
"""

import sys
import argparse
import numpy as np
import json
from pathlib import Path

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')

from flint import arb, ctx
import mpmath as mp


def build_arb_tau_components(c, N, prec):
    """
    Build τ = W02 - WR - Wp, returning each component separately for spectral analysis.

    Returns:
      (W02_matrix, WR_matrix, Wp_matrix, full_tau_matrix)
    """
    ctx.prec = prec

    from arb_ldlt_certify import arb_closed_forms, arb_J, arb_kappa, prime_powers_up_to
    from flint import arb_mat

    S, CC, XC, L = arb_closed_forms(N, c, prec)
    PI = arb.pi()
    sp2 = 16 * PI * PI
    l2 = L * L
    pref02 = 32 * L * (L / 4).sinh() ** 2
    kappa = arb_kappa(L)
    J = arb_J(L)

    pdata = prime_powers_up_to(c)
    weights = [arb(p).log() * (arb(q) ** arb('-0.5')) for (q, p) in pdata]
    positions = [arb(q).log() for (q, p) in pdata]
    M = len(weights)

    def S_signed(nn):
        return S[nn] if nn >= 0 else -S[-nn]

    DIM = 2 * N + 1
    W02_mat = arb_mat(DIM, DIM)
    WR_mat = arb_mat(DIM, DIM)
    Wp_mat = arb_mat(DIM, DIM)
    tau_mat = arb_mat(DIM, DIM)

    for i in range(DIM):
        n = i - N
        for j in range(i, DIM):
            m = j - N

            num = l2 - sp2 * m * n
            den = (l2 + sp2 * m * m) * (l2 + sp2 * n * n)
            W02_val = pref02 * num / den

            if n == m:
                WR_val = kappa + 2 * CC[abs(n)] + J - (2 / L) * XC[abs(n)]
            else:
                WR_val = (S_signed(m) - S_signed(n)) / (PI * (n - m))

            Wp_val = arb(0)
            for idx in range(M):
                y = positions[idx]
                if n == m:
                    q = 2 * (1 - y / L) * (2 * PI * n * y / L).cos()
                else:
                    q = ((2 * PI * m * y / L).sin() - (2 * PI * n * y / L).sin()) / (PI * (n - m))
                Wp_val += weights[idx] * q

            tau_val = W02_val - WR_val - Wp_val

            for mat, val in [(W02_mat, W02_val), (WR_mat, WR_val), (Wp_mat, Wp_val), (tau_mat, tau_val)]:
                mat[i, j] = val
                mat[j, i] = val

    return W02_mat, WR_mat, Wp_mat, tau_mat, DIM, M


def compute_spectral_properties(A_arb, DIM):
    """
    Convert Arb matrix to mpmath and compute eigenvalue spectrum.
    """
    ctx.prec = 300
    mp.mp.dps = 80

    Amp = mp.matrix(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            Amp[i, j] = mp.mpf(A_arb[i, j].str(200, radius=False))

    evals_mp = mp.eigsy(Amp, eigvals_only=True)
    evals = sorted([float(e) for e in evals_mp])
    return np.array(evals)


def analyze_universality(c, N, prec=300):
    """
    Extract properties relevant to universality class:
    1. Spectral norm of W02-WR (background)
    2. Spectral norm of Wp (primes)
    3. Number of negative eigenvalues of W02-WR
    4. Gap structure and condition number
    """
    print(f"\n" + "="*100)
    print(f"UNIVERSALITY CLASS ANALYSIS: c={c}, N={N}")
    print("="*100)

    W02, WR, Wp, tau, DIM, M = build_arb_tau_components(c, N, prec)

    # Compute W02-WR (background)
    mp.mp.dps = 80
    bg = np.zeros((DIM, DIM))
    for i in range(DIM):
        for j in range(DIM):
            val_str = (W02[i, j] - WR[i, j]).str(70, radius=False)
            bg[i, j] = float(mp.mpf(val_str))

    # Compute Wp
    wp_mat = np.zeros((DIM, DIM))
    for i in range(DIM):
        for j in range(DIM):
            val_str = Wp[i, j].str(70, radius=False)
            wp_mat[i, j] = float(mp.mpf(val_str))

    print(f"\nMatrix dimensions and rank:")
    print(f"  dim(τ) = {DIM}")
    print(f"  M = # prime powers ≤ {c} = {M}")
    print(f"  M/N = {M / N:.2f}  (prime density ratio)")

    # Spectral properties of background
    evals_bg = np.linalg.eigvalsh(bg)
    evals_bg_sorted = np.sort(evals_bg)
    n_neg_bg = np.sum(evals_bg_sorted < -1e-10)

    print(f"\nBackground matrix (W02 - WR) spectral properties:")
    print(f"  λ_min = {evals_bg_sorted[0]:+.6e}")
    print(f"  λ_max = {evals_bg_sorted[-1]:+.6e}")
    print(f"  ||·||_2 (spectral norm) = {np.max(np.abs(evals_bg_sorted)):+.6e}")
    print(f"  n_neg = {n_neg_bg}")
    print(f"  condition number κ(W02-WR) = {np.abs(evals_bg_sorted[-1]/evals_bg_sorted[0]) if evals_bg_sorted[0] != 0 else np.inf:.2e}")

    # Spectral properties of Wp
    evals_wp = np.linalg.eigvalsh(wp_mat)
    evals_wp_sorted = np.sort(evals_wp)

    print(f"\nPrime perturbation (Wp) spectral properties:")
    print(f"  λ_min(Wp) = {evals_wp_sorted[0]:+.6e}")
    print(f"  λ_max(Wp) = {evals_wp_sorted[-1]:+.6e}  (spectral norm)")
    print(f"  rank(Wp) ≈ {M} (number of prime powers, structural rank)")

    # Ratio of perturbation to background negativity
    neg_bulk = evals_bg_sorted[n_neg_bg - 1] if n_neg_bg > 0 else 0
    print(f"\nPerturbation-to-background ratio:")
    print(f"  |λ_min(W02-WR)| = {abs(evals_bg_sorted[0]):.6e}")
    print(f"  ||Wp||_2 = {evals_wp_sorted[-1]:.6e}")
    print(f"  ratio = ||Wp||_2 / |λ_min(BG)| = {evals_wp_sorted[-1] / abs(evals_bg_sorted[0]):.4f}")

    # Full matrix properties
    tau_np = np.zeros((DIM, DIM))
    for i in range(DIM):
        for j in range(DIM):
            val_str = tau[i, j].str(70, radius=False)
            tau_np[i, j] = float(mp.mpf(val_str))
    evals_tau = np.linalg.eigvalsh(tau_np)
    evals_tau_sorted = np.sort(evals_tau)

    print(f"\nFull matrix (τ = W02-WR-Wp) spectral properties:")
    print(f"  λ_min = {evals_tau_sorted[0]:+.6e}")
    print(f"  λ_max = {evals_tau_sorted[-1]:+.6e}")
    print(f"  n_neg = {np.sum(evals_tau_sorted < -1e-10)}")

    # Test hypothesis: scaling with c and N
    print(f"\nEmpirical scaling checks:")
    print(f"  log₁₀(λ_min) = {np.log10(max(abs(evals_tau_sorted[0]), 1e-20)):.3f}")
    print(f"  c^(-0.44N) ≈ 10^({-0.44 * N * np.log10(c):.2f})")
    print(f"  agreement: within 10x? {abs(np.log10(max(abs(evals_tau_sorted[0]), 1e-20))) - (-0.44 * N * np.log10(c)) < 1}")

    return {
        'c': c, 'N': N,
        'M': M,
        'bg_lambda_min': float(evals_bg_sorted[0]),
        'bg_lambda_max': float(evals_bg_sorted[-1]),
        'bg_n_neg': int(n_neg_bg),
        'wp_spectral_norm': float(evals_wp_sorted[-1]),
        'tau_lambda_min': float(evals_tau_sorted[0]),
        'ratio_wp_to_bg': float(evals_wp_sorted[-1] / abs(evals_bg_sorted[0])),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--c', type=int, default=30, help='cutoff parameter')
    p.add_argument('--N', type=int, default=8, help='number of frequency modes')
    p.add_argument('--prec', type=int, default=300, help='Arb precision in bits')
    args = p.parse_args()

    result = analyze_universality(args.c, args.N, prec=args.prec)

    # Save result
    output_json = f'universality_analysis_c{args.c}_N{args.N}.json'
    with open(output_json, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {output_json}")


if __name__ == '__main__':
    main()

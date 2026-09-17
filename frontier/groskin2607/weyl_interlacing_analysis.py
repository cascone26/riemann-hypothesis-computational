#!/usr/bin/env python3
"""
Weyl interlacing analysis for the parity-decomposed Weil form.

Key insight: τ = W02-WR - Wp decomposes into:
  τ_even = (W02-WR)_even - Wp_even
  τ_odd = (W02-WR)_odd - Wp_odd

Both sectors are PSD (verified), but the backgrounds are NOT:
  (W02-WR)_even has n_neg_even negative eigenvalues
  (W02-WR)_odd has n_neg_odd negative eigenvalues

Weyl's inequality bounds how much a rank-k perturbation can affect eigenvalues:
  λ_i(A + B) ≤ λ_i(A) + λ_1(B)  (if B is PSD)
  λ_i(A) + λ_n(B) ≤ λ_i(A+B)    (if B is any matrix)

Strategy: model Wp_even and Wp_odd as rank-limited perturbations and use interlacing
to bound the minimum eigenvalue of τ below by something that decays with (c, N).

The empirical law from margin_law.py is:
  λ_min ~ K * c^(-0.44*N)

We want to derive (or at least explain) this exponent analytically.
"""
import numpy as np
from scipy import linalg
import sys

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from cascade_generalizes_finiteN import build_tau_and_terms

def analyze_spectral_gaps(label, mat):
    """Analyze eigenvalue spectrum for gap structure."""
    evals = np.linalg.eigvalsh(mat)
    evals = np.sort(evals)[::-1]  # Sort descending

    gaps = np.diff(evals)
    max_gap_idx = np.argmax(gaps)
    max_gap = gaps[max_gap_idx]

    print(f"  {label}:")
    print(f"    evals (top 10): {evals[:10]}")
    print(f"    max gap: {max_gap:.6e} between evals[{max_gap_idx}] and evals[{max_gap_idx+1}]")

    return evals, max_gap

def rank_analysis_per_sector(c, N, W02_WR, Wp_terms, even_idx, odd_idx):
    """Analyze rank structure of Wp within each parity sector."""

    # Extract Wp perturbations per sector
    Wp_even_terms = []
    Wp_odd_terms = []

    for Wp in Wp_terms:
        Wp_even = Wp[np.ix_(even_idx, even_idx)]
        Wp_odd = Wp[np.ix_(odd_idx, odd_idx)]
        Wp_even_terms.append(Wp_even)
        Wp_odd_terms.append(Wp_odd)

    Wp_even_total = sum(Wp_even_terms)
    Wp_odd_total = sum(Wp_odd_terms)

    # SVD analysis
    _, s_even, _ = linalg.svd(Wp_even_total)
    _, s_odd, _ = linalg.svd(Wp_odd_total)

    def effective_rank(s, tol=1e-10):
        """Effective rank as number of singular values above threshold."""
        return np.sum(s > np.max(s) * tol)

    rank_even = effective_rank(s_even)
    rank_odd = effective_rank(s_odd)

    print(f"\nRank analysis (c={c}, N={N}):")
    print(f"  Wp_even: rank ≈ {rank_even} out of {len(even_idx)}, σ_1={s_even[0]:.6e}, σ_min={s_even[-1]:.6e}")
    print(f"  Wp_odd:  rank ≈ {rank_odd} out of {len(odd_idx)}, σ_1={s_odd[0]:.6e}, σ_min={s_odd[-1]:.6e}")

    return rank_even, rank_odd

def interlacing_bound(c, N, tau_even, tau_odd, W02WR_even, W02WR_odd, Wp_even, Wp_odd):
    """
    Apply Weyl interlacing to bound λ_min(τ_sector) in terms of structure.

    τ_sector = (W02-WR)_sector - Wp_sector

    By the Courant-Fischer minimax theorem:
    λ_min(A - B) >= λ_min(A) - λ_max(B)  (always holds)

    But for EIGENVALUE INTERLACING, if we know the spectrum of A,
    and B has effective rank r_eff, then only r_eff eigenvalues of A+B
    can deviate by more than the top singular value of B.

    More precisely (Weyl):
    λ_j(A+B) <= λ_j(A) + σ_1(B)
    """

    evals_even = np.linalg.eigvalsh(W02WR_even)
    evals_odd = np.linalg.eigvalsh(W02WR_odd)

    sv_even = np.linalg.svd(Wp_even, compute_uv=False)
    sv_odd = np.linalg.svd(Wp_odd, compute_uv=False)

    sigma_max_even = sv_even[0]
    sigma_max_odd = sv_odd[0]

    # Count how many evals(W02-WR) are below -sigma_max(Wp)
    # These MUST stay negative (can't be flipped by a rank-k perturbation)
    protected_even = np.sum(evals_even < -sigma_max_even)
    protected_odd = np.sum(evals_odd < -sigma_max_odd)

    print(f"\nInterlacing analysis (c={c}, N={N}):")
    print(f"  Even sector:")
    print(f"    (W02-WR)_even spectrum: min={evals_even[0]:.6e}, max={evals_even[-1]:.6e}")
    print(f"    σ_max(Wp_even) = {sigma_max_even:.6e}")
    print(f"    Protected negative: eigenvalues < {-sigma_max_even:.6e}: {protected_even} out of {len(evals_even)}")
    print(f"  Odd sector:")
    print(f"    (W02-WR)_odd spectrum: min={evals_odd[0]:.6e}, max={evals_odd[-1]:.6e}")
    print(f"    σ_max(Wp_odd) = {sigma_max_odd:.6e}")
    print(f"    Protected negative: eigenvalues < {-sigma_max_odd:.6e}: {protected_odd} out of {len(evals_odd)}")

    # Potential bound on λ_min(τ):
    # If τ is PD (both sectors PD), then λ_min(τ) >= min(λ_min(τ_even), λ_min(τ_odd))
    # But how small can λ_min(τ_sector) be?
    # By Courant-Fischer: λ_min(τ_sector) >= λ_{d_sector}(W02-WR)_sector - σ_max(Wp_sector)
    # where d_sector is the dimension of the sector

    # Actually, a better bound: the smallest eigenvalue of τ_sector
    # is at least the smallest "protected" eigenvalue (can't be flipped)
    # minus the perturbation magnitude. But if there are no protected ones,
    # it could be very small.

    bound_even = evals_even[protected_even] if protected_even > 0 else evals_even[-1] - sigma_max_even
    bound_odd = evals_odd[protected_odd] if protected_odd > 0 else evals_odd[-1] - sigma_max_odd

    print(f"\n  Potential lower bounds on λ_min(τ_sector):")
    print(f"    Even: {bound_even:.6e}")
    print(f"    Odd:  {bound_odd:.6e}")

    # Compute actual λ_min
    tau_even = W02WR_even - Wp_even
    tau_odd = W02WR_odd - Wp_odd
    lam_even = np.linalg.eigvalsh(tau_even)[0]
    lam_odd = np.linalg.eigvalsh(tau_odd)[0]

    print(f"\n  Actual λ_min(τ_sector):")
    print(f"    Even: {lam_even:.6e}")
    print(f"    Odd:  {lam_odd:.6e}")

    # The scaling law: empirically, λ_min ~ K * c^(-0.44*N)
    # Log analysis
    import math
    if lam_even > 1e-300:
        log_lam_even = math.log(abs(lam_even))
        print(f"\n  Even sector: log(|λ_min|) = {log_lam_even:.4f}")
        print(f"    vs log(c^(-0.44*N)) = {-0.44 * N * math.log(c):.4f}")

    return lam_even, lam_odd


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--c', type=int, default=30)
    parser.add_argument('--N', type=int, default=8)
    args = parser.parse_args()

    c, N = args.c, args.N
    dim = 2*N + 1

    print("="*100)
    print(f"Weyl Interlacing Analysis for the Parity-Decomposed Weil Form")
    print(f"c={c}, N={N}, dim={dim}")
    print("="*100)

    # Build matrices
    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)

    # Parity decomposition
    even_idx = [i for i in range(dim) if (i - N) % 2 == 0]
    odd_idx = [i for i in range(dim) if (i - N) % 2 == 1]

    W02WR_even = W02_WR[np.ix_(even_idx, even_idx)]
    W02WR_odd = W02_WR[np.ix_(odd_idx, odd_idx)]

    Wp_even = Wp_total[np.ix_(even_idx, even_idx)]
    Wp_odd = Wp_total[np.ix_(odd_idx, odd_idx)]

    # Spectral analysis
    print("\n--- SPECTRAL GAP ANALYSIS ---")
    evals_arch_even, gap_even = analyze_spectral_gaps(
        "(W02-WR)_even", W02WR_even)
    evals_arch_odd, gap_odd = analyze_spectral_gaps(
        "(W02-WR)_odd", W02WR_odd)

    # Rank analysis
    rank_even, rank_odd = rank_analysis_per_sector(c, N, W02_WR, Wp_terms, even_idx, odd_idx)

    # Interlacing analysis
    lam_even, lam_odd = interlacing_bound(
        c, N, W02WR_even - Wp_even, W02WR_odd - Wp_odd,
        W02WR_even, W02WR_odd, Wp_even, Wp_odd)

    print("\n" + "="*100)
    print("SUMMARY & OBSERVATIONS")
    print("="*100)
    print("""
Both even and odd sectors achieve PD status through a precision cancellation:
the negative eigenvalue spectrum of (W02-WR)_sector is exactly matched and cancelled
by the Wp_sector perturbation.

This suggests:
1. The perturbation Wp is NOT rank-limited in a way that simple interlacing bounds
   (which assume rank-k perturbations) can directly apply.

2. Instead, the cancellation is more subtle: each prime power q<=c contributes a
   rank-2 matrix (from sin/cos of 2*pi*n*log(q)/L), and these stack up to collectively
   cancel the arch negativity.

3. The 0.44 exponent likely emerges from:
   - The number of prime powers M ~ c/log(c) (grows with c)
   - The dimension 2N+1 (grows with N)
   - A precise cancellation law that becomes LESS precise (higher margin) as
     one or the other grows.

Next steps to refine the bound:
(A) Understand the per-prime-power rank-2 contribution structure more explicitly.
(B) Try a tighter interlacing bound that accounts for MULTIPLE rank-2 blocks,
    not just the top singular value.
(C) Compare with the empirical cascading Hellmann-Feynman decomposition from S2624
    to see if it reveals the actual cancellation structure.
""")

if __name__ == '__main__':
    main()

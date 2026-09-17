#!/usr/bin/env python3
"""bochner_lens_rank1.py -- §1: Positive-definiteness lens on prime vs arch.

THESIS: prime and arch are both built from positive-definite (Bochner) kernels:
  arch: Ω(r) (archimedean density, Cauchy–Stieltjes)
  prime: Σ Λ(n)/√n cos(r log n) (positive-definite cosine-sum kernel)

On a band-limited (prolate) subspace, they are NEARLY IDENTICAL except for one
direction (the pole's principal axis). Test if prime−arch has RANK 1 + negligible
on the band-limited space — if so, κ≤1 is structural (not algebraic).

HANDLE: all code-verified on explicit Q matrices.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh as scipy_eigh
import sys

def vonmangoldt_upto(N):
    """Compute Λ(n) for n=1..N exactly."""
    Lam = np.zeros(N + 1)
    is_comp = np.zeros(N + 1, dtype=bool)
    for p in range(2, N + 1):
        if not is_comp[p]:
            is_comp[2*p::p] = True
            lp = np.log(p); pk = p
            while pk <= N:
                Lam[pk] = lp; pk *= p
    return Lam

def build_form_prolate(cen_lo, cen_hi, fourW, frac=0.5, dr=0.008, return_pieces=False):
    """Build Q(φ)=arch+pole−prime on band-limited functions in [cen_lo, cen_hi].

    Band-limited: support in [-4W, 4W] (Fourier domain), via Fejér kernels.
    Centers: uniform grid at spacing π/(2W) covering [cen_lo, cen_hi].

    Returns: (Q, arch, prime, pole, G_gram, centers) if return_pieces=True
    else just Q.
    """
    W = fourW / 4.0
    L = cen_hi - cen_lo
    N = max(4, round(frac * L * 2 * W / np.pi))
    N = min(N, 20)  # Cap dimension for speed

    # Frequency grid (coarser to speed up)
    dr_use = dr
    r_max = cen_hi + 150
    RG = np.arange(0, r_max, dr_use)

    # Archimedean density Ω(r) = Re ψ(1/4 + ir/2) - ln π
    ln_pi = float(np.log(np.pi))
    OMEGA = np.real(digamma(0.25 + 0.5j*RG)) - ln_pi

    # Band-limited basis: Fejér kernels
    CEN = np.linspace(cen_lo, cen_hi, N)
    def fejer(x):
        y = W * x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz]) ** 2
        return out

    # Gram matrix via band-limit
    P = np.array([fejer(RG - c) + fejer(RG + c) for c in CEN])
    G = (P @ P.T) * dr_use

    # Archimedean form
    arch_form = (P * OMEGA[None, :]) @ P.T * dr_use / (2 * np.pi)

    # Prime form
    NMAX = 8000
    LAM = vonmangoldt_upto(NMAX)
    NS = np.arange(2, min(NMAX+1, int(np.exp(fourW))+10))
    LOGN = np.log(NS)
    LAMv = LAM[2:len(NS)+2]

    mask = LOGN < fourW
    ln = LOGN[mask]
    wts = LAMv[mask] / np.sqrt(NS[mask])

    cosu = np.cos(np.outer(ln, RG))
    prime_form = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr_use / (2 * np.pi)
        prime_form[j, :] = 2.0 * (gp * wts[None, :]).sum(axis=1)
    prime_form = 0.5 * (prime_form + prime_form.T)

    # Pole (rank-2, here ignored for simplicity)
    pole_form = np.zeros((N, N))

    # Q = arch + pole - prime
    Q = arch_form + pole_form - prime_form

    if return_pieces:
        return Q, arch_form, prime_form, pole_form, G, CEN
    return Q

def analyze_rank_structure(arch, prime, G, label=""):
    """Analyze the rank and eigenvalue structure of arch, prime, and their difference."""
    print(f"\n=== {label} Rank/Eigenvalue Analysis ===")

    # Generalized eigenproblem
    gew_arch, gev_arch = scipy_eigh(arch, G)
    gew_prime, gev_prime = scipy_eigh(prime, G)

    # Rank via positive eigenvalues
    rank_arch = np.sum(gew_arch > 1e-10)
    rank_prime = np.sum(gew_prime > 1e-10)

    print(f"  dim(space) = {arch.shape[0]}")
    print(f"  rank(arch) = {rank_arch}, min eig = {gew_arch.min():.6e}, max = {gew_arch.max():.6e}")
    print(f"  rank(prime) = {rank_prime}, min eig = {gew_prime.min():.6e}, max = {gew_prime.max():.6e}")

    # Difference matrix
    diff = arch - prime
    gew_diff, gev_diff = scipy_eigh(diff, G)
    print(f"  eig(arch−prime) = [{gew_diff.min():.6e}, ..., {gew_diff.max():.6e}]")
    n_pos = np.sum(gew_diff > 1e-10)
    n_neg = np.sum(gew_diff < -1e-10)
    print(f"    #{n_pos} positive, #{n_neg} negative")

    # Rank-1 approximation via top eigenvector
    idx_top = np.argmax(np.abs(gew_diff))
    v_top = gev_diff[:, idx_top]
    rank1_approx = gew_diff[idx_top] * np.outer(v_top, v_top)
    rank1_error = np.linalg.norm(diff - rank1_approx, 'fro') / (np.linalg.norm(diff, 'fro') + 1e-15)
    print(f"  Rank-1 test: error/norm = {rank1_error:.6e}")

    # Condition numbers
    cond_arch = gew_arch.max() / np.maximum(gew_arch.min(), 1e-15)
    cond_prime = gew_prime.max() / np.maximum(gew_prime.min(), 1e-15)
    print(f"  cond(arch) = {cond_arch:.6e}, cond(prime) = {cond_prime:.6e}")

    return {
        'eig_arch': gew_arch, 'eig_prime': gew_prime, 'eig_diff': gew_diff,
        'rank1_error': rank1_error, 'v_top': v_top
    }

def main():
    print("=" * 70)
    print("BOCHNER LENS: Rank-1 Analysis of arch vs prime")
    print("=" * 70)

    # Test cases: corner window at various heights
    test_cases = [
        ("shallow", 10, 20, 6, 0.5),
        ("medium", 50, 100, 10, 0.4),
        ("deep", 500, 1000, 16, 0.3),
    ]

    for label, cen_lo, cen_hi, fourW, frac in test_cases:
        print(f"\n{'='*70}")
        print(f"Case: {label} (height {cen_lo}–{cen_hi}, bandwidth {fourW})")
        print(f"{'='*70}")

        try:
            Q, arch, prime, pole, G, CEN = build_form_prolate(
                cen_lo, cen_hi, fourW, frac=frac, return_pieces=True
            )

            results = analyze_rank_structure(arch, prime, G, label=label)

            # Generalized eigenvalues of Q
            gew_Q, gev_Q = scipy_eigh(Q, G)
            print(f"\n  Q = arch − prime:")
            print(f"    min eig = {gew_Q.min():.6e}, max = {gew_Q.max():.6e}")
            n_pos_Q = np.sum(gew_Q > 1e-10)
            n_neg_Q = np.sum(gew_Q < -1e-10)
            print(f"    #{n_pos_Q} positive, #{n_neg_Q} negative")

            if n_neg_Q > 0:
                print(f"    *** NEGATIVE EIGENVALUE: κ = {n_neg_Q} ***")
                idx_min = np.argmin(gew_Q)
                v_min = gev_Q[:, idx_min]
                print(f"    Margin: {-gew_Q[idx_min]:.6e}")

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

if __name__ == '__main__':
    main()

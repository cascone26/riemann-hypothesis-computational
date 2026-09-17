#!/usr/bin/env python3
"""variational_extremal.py -- Rayleigh quotient via explicit arch/prime construction.

LENS: Build arch and prime explicitly in the grid-basis (like kappa_weyl_count.py),
then study the generalized eigenvalue problem:

  prime x = λ arch x  ⟹  κ = #{λ > 1}

The variational interpretation:

  λ_max = max_g (g^T prime g) / (g^T arch g)

If arch > 0 and prime > 0 (both PSD), then λ > 0.
For RH, we need λ ≤ 1 for all but one direction (the pole).

HYPOTHESIS: Can we unconditionally bound λ_2 ≤ 1 from structural properties of arch?
That would give κ ≤ 1 without using the pole-vs-negative alignment.
"""

import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

def von_mangoldt(N):
    """Compute Λ(n) for n=1..N."""
    lam = np.zeros(N + 1)
    is_prime = np.ones(N + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(np.sqrt(N)) + 1):
        if is_prime[p]:
            is_prime[p*p::p] = False
    for p in range(2, N + 1):
        if is_prime[p]:
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

def build_arch_prime_explicit(cen_lo, cen_hi, fourW, frac=0.8, dr=0.004):
    """
    Explicitly build arch and prime matrices (like §965's approach).

    Returns:
      arch: Archimedean quadratic form (should be PSD)
      prime: Prime quadratic form (should be PSD)
      G: The Gram matrix (for projecting to the kernel)
    """
    W = fourW / 4.0
    L = cen_hi - cen_lo
    N = max(4, round(frac * L * 2 * W / np.pi))

    # r-grid for quadrature
    RG = np.arange(-(cen_hi + 120), (cen_hi + 120), dr)

    # Archimedean density ψ = Re(ψ(1/4 + ir/2)) - log π
    PSI_val = np.real(digamma(0.25 + 0.5j * RG)) - np.log(np.pi)

    # Fejér centers
    CEN = np.linspace(cen_lo, cen_hi, N)

    # Fejér window basis
    def fejér(x):
        y = W * x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz]) ** 2
        return out

    # Fejér basis: P[j] = fejér(r - c_j) + fejér(r + c_j)
    P = np.array([fejér(RG - c) + fejér(RG + c) for c in CEN])

    # Gram matrix: G = ∫ P^T P
    G = (P @ P.T) * dr

    # Arch = ∫ ψ(r) P^T P
    arch = (P * PSI_val[None, :]) @ P.T * dr / (2 * np.pi)

    # Prime: ∫ (Σ Λ(n) cos(log(n) r)) P
    NMAX = max(1000, int(np.exp(fourW)))
    LAM = von_mangoldt(NMAX)
    NS = np.arange(2, min(NMAX + 1, 10000))
    LOGN = np.log(NS)
    LAMv = LAM[NS]

    mask = LOGN < fourW
    logn = LOGN[mask]
    wts = LAMv[mask] / np.sqrt(NS[mask])

    cosu = np.cos(np.outer(logn, RG))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2 * np.pi)
        prime[j, :] = 2.0 * (gp * wts[None, :]).sum(axis=1)

    # Symmetrize
    prime = 0.5 * (prime + prime.T)

    return arch, prime, G, N


def generalized_spectrum_analysis(arch, prime, G, verbose=True):
    """
    Analyze the generalized eigenvalue problem: prime x = λ arch x.

    Uses the Gram matrix to project to the space orthogonal to the null.
    """
    # Project to orthogonal complement of kernel (using Gram matrix)
    gw, gv = np.linalg.eigh(G)
    keep = gw > gw.max() * 1e-12
    B = gv[:, keep] / np.sqrt(gw[keep])

    # Compress to the non-kernel space
    arch_compressed = B.T @ arch @ B
    prime_compressed = B.T @ prime @ B

    # Check PSD
    arch_eigs = np.linalg.eigvalsh(arch_compressed)
    prime_eigs = np.linalg.eigvalsh(prime_compressed)

    if verbose:
        print(f"\nCompressed operators ({len(keep)} dims):")
        print(f"  arch:  neg={np.sum(arch_eigs < -1e-13)}, min={arch_eigs[0]:.3e}, max={arch_eigs[-1]:.3e}")
        print(f"  prime: neg={np.sum(prime_eigs < -1e-13)}, min={prime_eigs[0]:.3e}, max={prime_eigs[-1]:.3e}")

    # Solve generalized eigenvalue problem
    try:
        lam, X = eigh(prime_compressed, arch_compressed)
    except np.linalg.LinAlgError:
        if verbose:
            print("⚠ Generalized eigenvalue solve failed")
        return None

    kappa = np.sum(lam > 1.0 + 1e-10)

    if verbose:
        print(f"\nGENERALIZED SPECTRUM (prime x = λ arch x):")
        print(f"  κ = #{'{λ > 1}'} = {kappa}")
        print(f"  λ_max = {lam[-1]:.6e}  (λ_max - 1 = {lam[-1]-1:.6e})")
        if len(lam) >= 2:
            print(f"  λ_2   = {lam[-2]:.6e}  (1 - λ_2 = {1-lam[-2]:.6e})")
            print(f"  λ_2 - λ_3 = {lam[-2] - lam[-3]:.6e}")
        print(f"  min λ = {lam[0]:.6e}")

        print(f"\nSpectrum tail:")
        for i in range(max(0, len(lam)-8), len(lam)):
            marker = " ← λ_top" if i == len(lam)-1 else (" ← λ_2nd" if i == len(lam)-2 else "")
            print(f"  λ[{i}] = {lam[i]:12.6e}{marker}")

    return {'lambda': lam, 'X': X, 'kappa': kappa, 'arch_eigs': arch_eigs, 'prime_eigs': prime_eigs}


def variational_scan_windows(window_list=[(14, 44), (30, 60), (45, 75)], fourW=6.0, frac=0.8):
    """Scan different windows and track κ, λ_2."""

    print("\n" + "="*90)
    print("VARIATIONAL SCAN: κ and λ_2 across windows at fixed 4W={:.1f}".format(fourW))
    print("="*90)

    print(f"{'Window':>12} {'N':>4} | {'λ_max':>12} {'λ_max-1':>12} {'λ_2':>12} {'1-λ_2':>12} {'κ':>3}")
    print("-"*95)

    for a, b in window_list:
        arch, prime, G, N = build_arch_prime_explicit(a, b, fourW, frac=frac)
        res = generalized_spectrum_analysis(arch, prime, G, verbose=False)
        if res is not None:
            lam = res['lambda']
            lam_2 = lam[-2] if len(lam) >= 2 else np.nan
            print(f"[{a:3d},{b:3d}]  {N:>4} | {lam[-1]:12.6e} {lam[-1]-1:12.6e} "
                  f"{lam_2:12.6e} {1-lam_2:12.6e} {res['kappa']:3d}")


def bandwidth_resolution_threshold(a_lo=14, a_hi=44):
    """
    Key question: at what 4W does κ jump from 0 to 1?
    This relates to the explicit-formula prime-zero resolution threshold.
    """
    print("\n" + "="*90)
    print("RESOLUTION THRESHOLD: 4W needed for κ to jump 0→1 at window [{},{}]".format(a_lo, a_hi))
    print("="*90)

    print(f"{'4W':>6} {'N':>4} | {'κ':>3} {'λ_max':>12} {'λ_2':>12} {'λ_3':>12}")
    print("-"*70)

    for fourW in [3.0, 4.0, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]:
        arch, prime, G, N = build_arch_prime_explicit(a_lo, a_hi, fourW, frac=0.8)
        res = generalized_spectrum_analysis(arch, prime, G, verbose=False)
        if res is not None:
            lam = res['lambda']
            lam_2 = lam[-2] if len(lam) >= 2 else np.nan
            lam_3 = lam[-3] if len(lam) >= 3 else np.nan
            print(f"{fourW:6.1f} {N:4d} | {res['kappa']:3d} {lam[-1]:12.6e} "
                  f"{lam_2:12.6e} {lam_3:12.6e}")


def second_eigenvalue_structure_detailed(a_lo=14, a_hi=44, fourW=6.0):
    """Deep dive into λ_2: its localization, alignment, structure."""

    print("\n" + "="*90)
    print(f"DETAILED λ_2 ANALYSIS for window [{a_lo},{a_hi}], 4W={fourW}")
    print("="*90)

    arch, prime, G, N = build_arch_prime_explicit(a_lo, a_hi, fourW, frac=0.8)

    # Project to kernel complement
    gw, gv = np.linalg.eigh(G)
    keep = gw > gw.max() * 1e-12
    B = gv[:, keep] / np.sqrt(gw[keep])

    arch_c = B.T @ arch @ B
    prime_c = B.T @ prime @ B

    # Generalized spectrum
    from scipy.linalg import eigh as sp_eigh
    lam, X = sp_eigh(prime_c, arch_c)
    kappa = np.sum(lam > 1.0 + 1e-10)

    print(f"κ = {kappa}, dim (after compression) = {len(lam)}")
    print(f"\nTop eigenvalues:")
    for i in range(min(5, len(lam))):
        idx = len(lam) - 1 - i
        print(f"  λ[{idx}] = λ_{{top-{i}}} = {lam[idx]:.6e}")

    # Analyze λ_2: is it close to 1?
    if len(lam) >= 2:
        lam_2 = lam[-2]
        x_2nd = X[:, -2]

        gap_to_top = lam[-1] - lam_2
        gap_from_below = lam_2 - lam[-3] if len(lam) >= 3 else np.nan

        print(f"\nλ_2 analysis:")
        print(f"  λ_2 = {lam_2:.6e}")
        print(f"  1 - λ_2 = {1 - lam_2:.6e}  (margin from 1)")
        print(f"  λ_max - λ_2 = {gap_to_top:.6e}  (gap to max)")
        if not np.isnan(gap_from_below):
            print(f"  λ_2 - λ_3 = {gap_from_below:.6e}  (gap from below)")

        # Rayleigh quotient check
        R_2 = (x_2nd @ prime_c @ x_2nd) / (x_2nd @ arch_c @ x_2nd)
        print(f"  R(x_2nd) verification = {R_2:.6e}  (should = λ_2 = {lam_2:.6e})")

    return {'lambda': lam, 'X': X, 'arch': arch, 'prime': prime}


if __name__ == "__main__":
    print("\n" + "="*90)
    print("VARIATIONAL EXTREMAL ANALYSIS — Rayleigh quotient of (prime, arch)")
    print("="*90)

    # Main analysis
    second_eigenvalue_structure_detailed(a_lo=14, a_hi=44, fourW=6.0)

    # Scan different windows
    variational_scan_windows(fourW=6.0)

    # Resolution threshold
    bandwidth_resolution_threshold()

    print("\n✓ Variational analysis complete.")

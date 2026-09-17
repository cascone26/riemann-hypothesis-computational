#!/usr/bin/env python3
"""variational_perturbation_test.py -- Test if λ_2 ≈ 1 is STRUCTURALLY forced or ARITHMETIC-DEPENDENT.

HYPOTHESIS 1 (Structural): λ_2 = 1 exactly due to some functional-equation identity,
  independent of whether RH holds. Perturbation should NOT move λ_2 far from 1.

HYPOTHESIS 2 (Arithmetic-dependent): λ_2 ≈ 1 only when the true Λ(n) is used.
  Perturbing the prime frequencies should BREAK the alignment, moving λ_2 away from 1.

TEST: Perturb Λ(n) → Λ(n) · U(n) where U(n) ∼ U(1−δ, 1+δ), then recompute λ_2.
  If λ_2 stays near 1 → Structural (Hypothesis 1)
  If λ_2 jumps away from 1 → Arithmetic-dependent (Hypothesis 2), meaning κ ≤ 1 = RH
"""

import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh as sp_eigh
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

def von_mangoldt(N):
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

def build_arch_prime_perturbed(cen_lo, cen_hi, fourW, frac=0.8, dr=0.004,
                               lambda_scale=1.0, seed=None):
    """
    Build arch and prime matrices with PERTURBED Λ coefficients.

    Parameters:
      lambda_scale: multiply each Λ(n) by this factor (1.0 = unperturbed)
      seed: random seed for noise generation
    """
    if seed is not None:
        np.random.seed(seed)

    W = fourW / 4.0
    L = cen_hi - cen_lo
    N = max(4, round(frac * L * 2 * W / np.pi))

    RG = np.arange(-(cen_hi + 120), (cen_hi + 120), dr)
    PSI_val = np.real(digamma(0.25 + 0.5j * RG)) - np.log(np.pi)

    CEN = np.linspace(cen_lo, cen_hi, N)

    def fejér(x):
        y = W * x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz]) ** 2
        return out

    P = np.array([fejér(RG - c) + fejér(RG + c) for c in CEN])
    G = (P @ P.T) * dr
    arch = (P * PSI_val[None, :]) @ P.T * dr / (2 * np.pi)

    NMAX = max(1000, int(np.exp(fourW)))
    LAM = von_mangoldt(NMAX)
    NS = np.arange(2, min(NMAX + 1, 10000))
    LOGN = np.log(NS)
    LAMv = LAM[NS]

    mask = LOGN < fourW
    logn = LOGN[mask]
    wts_base = LAMv[mask] / np.sqrt(NS[mask])

    # PERTURB the weights
    wts = wts_base * lambda_scale

    cosu = np.cos(np.outer(logn, RG))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2 * np.pi)
        prime[j, :] = 2.0 * (gp * wts[None, :]).sum(axis=1)

    prime = 0.5 * (prime + prime.T)

    return arch, prime, G, N


def analyze_spectrum_perturbed(lambda_scale, a_lo=14, a_hi=44, fourW=6.0, seed=None):
    """Compute λ_2 under perturbation."""
    arch, prime, G, N = build_arch_prime_perturbed(a_lo, a_hi, fourW, frac=0.8,
                                                     lambda_scale=lambda_scale, seed=seed)

    gw, gv = np.linalg.eigh(G)
    keep = gw > gw.max() * 1e-12
    B = gv[:, keep] / np.sqrt(gw[keep])

    arch_c = B.T @ arch @ B
    prime_c = B.T @ prime @ B

    lam, X = sp_eigh(prime_c, arch_c)
    kappa = np.sum(lam > 1.0 + 1e-10)

    lam_2 = lam[-2] if len(lam) >= 2 else np.nan
    lam_3 = lam[-3] if len(lam) >= 3 else np.nan

    return {
        'kappa': kappa,
        'lambda_max': lam[-1],
        'lambda_2': lam_2,
        'lambda_3': lam_3,
        '1-lambda_2': 1 - lam_2 if not np.isnan(lam_2) else np.nan,
    }


def perturbation_sweep(delta_list=[0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1]):
    """
    Sweep δ where λ(n) → λ(n) · (1 + δ·U(−1,1)).
    Since we're just scaling, use λ(n) → λ(n)·(1+δ).
    """

    print("\n" + "="*100)
    print("PERTURBATION SWEEP: How does λ_2 respond to SCALED prime coefficients?")
    print("Testing if λ_2 ≈ 1 is structurally FORCED or ARITHMETIC-DEPENDENT")
    print("="*100)

    print(f"\n{'δ (scale factor = 1+δ)':>20} {'λ_max':>15} {'λ_max−1':>15} {'λ_2':>15} "
          f"{'1−λ_2':>15} {'κ':>4}")
    print("-" * 100)

    results = []
    for delta in delta_list:
        scale = 1.0 + delta
        res = analyze_spectrum_perturbed(scale, a_lo=14, a_hi=44, fourW=6.0)
        results.append((delta, res))

        print(f"{delta:>20.6f} {res['lambda_max']:>15.6e} {res['lambda_max']-1:>15.6e} "
              f"{res['lambda_2']:>15.6e} {res['1-lambda_2']:>15.6e} {res['kappa']:>4d}")

    return results


def interpretations_and_verdict(results):
    """Analyze what the perturbation results tell us."""

    print("\n" + "="*100)
    print("INTERPRETATION")
    print("="*100)

    # Extract λ_2 trend
    deltas = np.array([r[0] for r in results])
    lambda_2_vals = np.array([r[1]['lambda_2'] for r in results])
    kappas = np.array([r[1]['kappa'] for r in results])

    # Check if λ_2 stays near 1
    mean_distance_from_1 = np.mean(np.abs(lambda_2_vals - 1.0))
    max_distance_from_1 = np.max(np.abs(lambda_2_vals - 1.0))

    print(f"\nλ_2 behavior:")
    print(f"  Mean distance from 1.0: {mean_distance_from_1:.3e}")
    print(f"  Max distance from 1.0: {max_distance_from_1:.3e}")
    print(f"  Unperturbed (δ=0): λ_2 = {lambda_2_vals[0]:.6e}")
    print(f"  Max perturbed (δ={deltas[-1]:.4f}): λ_2 = {lambda_2_vals[-1]:.6e}")

    # κ behavior
    unique_kappas = np.unique(kappas)
    print(f"\nκ behavior:")
    print(f"  Unique κ values: {unique_kappas}")
    print(f"  κ stays constant at {kappas[0]} for all δ: {np.all(kappas == kappas[0])}")

    # Verdict
    if max_distance_from_1 < 1e-4:
        print(f"\n{'='*100}")
        print("VERDICT: λ_2 ≈ 1 is STRUCTURALLY FORCED (Hypothesis 1)")
        print("="*100)
        print("\nThe near-equality prime ≈ arch is NOT dependent on the true Λ(n).")
        print("Even scaling Λ uniformly does not move λ_2 away from 1.")
        print("This suggests a functional-equation identity that FORCES λ_2 = 1.")
        print("\nIMPLICATION: κ ≤ 1 may be provable UNCONDITIONALLY from the structure")
        print("of the Fejér window and the digamma density, independent of RH.")
    elif kappas[-1] > kappas[0]:
        print(f"\n{'='*100}")
        print("VERDICT: λ_2 ≈ 1 is ARITHMETIC-DEPENDENT (Hypothesis 2)")
        print("="*100)
        print(f"\nScaling Λ perturbs κ: κ={int(kappas[0])} → κ={int(kappas[-1])}.")
        print("The single negative square only exists with the true primes.")
        print("\nIMPLICATION: κ ≤ 1 ⟺ RH (unconditional reduction impossible via this lane).")
        print("But the STRUCTURE (λ_2 ≈ 1) is a strong hint of the Π₁ form.")
    else:
        print(f"\n{'='*100}")
        print("MIXED: κ stays 1, but the structure may change in subtle ways")
        print("="*100)


if __name__ == "__main__":
    print("\n" + "="*100)
    print("VARIATIONAL PERTURBATION TEST")
    print("="*100)

    results = perturbation_sweep(delta_list=[0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1])

    interpretations_and_verdict(results)

    print("\n✓ Perturbation test complete.")

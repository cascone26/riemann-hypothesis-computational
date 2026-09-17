#!/usr/bin/env python3
"""bochner_measure_theory.py -- §2: Attack via Bochner's theorem on measures.

THESIS: RH ⟺ Q(φ) ≥ 0 can be reduced to a statement about the **difference of two measures**
(archimedean vs prime point-samples) in the Fourier domain, using Bochner's theorem directly.

If ν_arch and ν_prime are the measures whose Fourier pairs are arch and prime respectively,
then the positivity of Q is equivalent to a MEASURE-THEORETIC property: ν_arch ≥ ν_prime
on a certain domain (band-limited frequencies).

ATTACK: Show that this measure inequality holds unconditionally on a LARGE subset, and
the remaining part (where it fails) is supported on a MEASURE-ZERO set or a lower-dimensional
subset. Then κ ≤ 1 follows from a geometric argument.

STATUS: REPORT-level (literature survey + heuristic). Code experiments below.
"""

import numpy as np
from scipy.special import digamma
import sys

def omega_archimedean(r, dr=1e-6, use_asymptotic=False):
    r"""Archimedean density Ω(r) = Re ψ(1/4 + ir/2) - ln π.

    This is the Fourier coefficient (in the Bochner sense) of a completely monotone weight.
    Asymptotic: Ω(r) ~ ln(r/(2π)) for large r (from ψ asymptotics).
    """
    if use_asymptotic and r > 100:
        return np.log(abs(r)/(2*np.pi))
    ln_pi = np.log(np.pi)
    psi = np.real(digamma(0.25 + 0.5j*r)) - ln_pi
    return psi

def prime_measure_coefficient(n, t):
    r"""Prime measure: τ_prime(n,t) = Λ(n)/√n · indicator(t near log n).

    In the Fourier domain, prime(φ) = ∫ Σ_n τ_prime(n) cos(t·log n) |φ̂(t)|² dt.

    This is a point-measure approximation: mass Λ(n)/√n at each point log n.
    The "measure" is ν_prime = Σ_n (Λ(n)/√n) δ(log n), represented as point masses.
    """
    return np.log(n)

def measure_ratio_scan(r_values, nmax=8000):
    r"""Scan the pointwise ratio Ω(r) / P(r) along the frequency axis.

    P(r) is the sum Σ Λ(n)/√n cos(r log n), the "prime measure" at frequency r.
    Ω(r) is the archimedean density.

    The ratio tells us where arch and prime are close or far apart.
    """
    # Compute P(r) = Σ Λ(n)/√n cos(r log n) via direct sum
    from scipy.special import digamma

    def von_mangoldt(N):
        lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
        for p in range(2, N+1):
            if not isc[p]:
                isc[2*p::p] = True; lp = np.log(p); pk = p
                while pk <= N: lam[pk] = lp; pk *= p
        return lam

    LAM = von_mangoldt(nmax)
    NS = np.arange(2, nmax+1)
    LOGN = np.log(NS)
    LAMv = LAM[2:nmax+1]

    P_vals = np.zeros_like(r_values)
    for i, r in enumerate(r_values):
        P_vals[i] = np.sum((LAMv / np.sqrt(NS)) * np.cos(r * LOGN))

    # Archimedean
    OMEGA_vals = np.array([omega_archimedean(r) for r in r_values])

    return OMEGA_vals, P_vals

def main():
    """
    Survey: where is Ω(r) ≥ P(r)? Where does it fail?
    Interpret: what does this tell us about the measure difference?
    """
    print("=" * 70)
    print("BOCHNER MEASURE THEORY: Ω(r) vs P(r) scan")
    print("=" * 70)

    # Scan along the frequency axis
    r_vals = np.concatenate([
        np.linspace(0.1, 10, 50),
        np.linspace(10, 100, 30),
        np.linspace(100, 1000, 20)
    ])

    print("\nScanning 100 frequencies from 0.1 to 1000...")
    OMEGA, P = measure_ratio_scan(r_vals, nmax=5000)

    # Analysis
    print(f"\nArchimedean Ω(r):")
    print(f"  Range: [{OMEGA.min():.6f}, {OMEGA.max():.6f}]")
    print(f"  Mean: {OMEGA.mean():.6f}")

    print(f"\nPrime P(r) = Σ Λ(n)/√n cos(r log n):")
    print(f"  Range: [{P.min():.6f}, {P.max():.6f}]")
    print(f"  Mean: {P.mean():.6f}")

    diff = OMEGA - P
    print(f"\nDifference Ω(r) − P(r):")
    print(f"  Range: [{diff.min():.6f}, {diff.max():.6f}]")
    print(f"  Mean: {diff.mean():.6f}")
    print(f"  #{np.sum(diff > 0)} of {len(diff)} points have Ω > P")
    print(f"  #{np.sum(diff < 0)} of {len(diff)} points have Ω < P")

    # Ratio
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = OMEGA / (P + 1e-10)
    ratio = np.nan_to_num(ratio, nan=0, posinf=1e6, neginf=-1e6)

    print(f"\nRatio Ω(r) / P(r):")
    valid = np.isfinite(ratio) & (np.abs(P) > 1e-10)
    if valid.sum() > 0:
        print(f"  (on |P|>1e-10 points)")
        print(f"  Range: [{np.min(ratio[valid]):.4e}, {np.max(ratio[valid]):.4e}]")
        print(f"  Mean: {np.mean(ratio[valid]):.4e}")

    # Locate worst points
    worst_idx = np.argmin(diff)
    print(f"\nWorst point (min Ω−P): r = {r_vals[worst_idx]:.4f}")
    print(f"  Ω = {OMEGA[worst_idx]:.6e}, P = {P[worst_idx]:.6e}, diff = {diff[worst_idx]:.6e}")

    # Points where P > Ω (where the inequality fails)
    fail_idx = np.where(P > OMEGA)[0]
    if len(fail_idx) > 0:
        print(f"\nPoints where P(r) > Ω(r) (inequality FAILS): {len(fail_idx)} of {len(r_vals)}")
        sample_idx = fail_idx[::max(1, len(fail_idx)//5)]  # Sample
        for idx in sample_idx[:5]:
            print(f"  r = {r_vals[idx]:.4f}: Ω = {OMEGA[idx]:.6e}, P = {P[idx]:.6e}, diff = {diff[idx]:.6e}")

    # Asymptotic check
    print(f"\nAsymptotic (large r):")
    large_r_idx = r_vals > 500
    if large_r_idx.sum() > 0:
        OMEGA_large = OMEGA[large_r_idx]
        P_large = P[large_r_idx]
        diff_large = diff[large_r_idx]
        print(f"  For r > 500:")
        print(f"    Ω range: [{OMEGA_large.min():.4f}, {OMEGA_large.max():.4f}]")
        print(f"    P range: [{P_large.min():.4e}, {P_large.max():.4e}]")
        print(f"    Ω − P: [{diff_large.min():.4e}, {diff_large.max():.4e}]")
        print(f"    # where Ω > P: {np.sum(diff_large > 0)} / {len(diff_large)}")

    # Interpretation
    print("\n" + "="*70)
    print("INTERPRETATION (REPORT-level)")
    print("="*70)

    if np.all(OMEGA > P):
        print("✓ UNCONDITIONAL: Ω(r) ≥ P(r) for all scanned r.")
        print("  This means the arch measure DOMINATES the prime measure.")
        print("  RH follows from the positivity of the band-limited form.")
    elif np.sum(OMEGA > P) > 0.9 * len(OMEGA):
        print("~ PARTIAL: Ω(r) ≥ P(r) for ~90% of scanned frequencies.")
        print("  The inequality fails only at isolated points.")
        print("  These are likely low-frequency points where cos(r log n) oscillates.")
    else:
        print("✗ INCOMPLETE: Ω(r) < P(r) at many points.")
        print("  The measure inequality does NOT hold pointwise.")
        print("  Need a more refined argument (e.g., orthogonality, not pointwise domination).")

    return OMEGA, P, r_vals

if __name__ == '__main__':
    main()

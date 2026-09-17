#!/usr/bin/env python3
"""
gaussian_explicit_formula_check.py — verify-before-ship cross-check for the Gaussian family.

The Gaussian scan showed tau(w)=A+B-P underflows for w<~2.3 because A+B-P IS the Weil explicit
formula = Σ_γ h(γ) with h(r)=exp(-r^2/w^2), and the first zero γ1=14.13 gives exp(-γ1^2/w^2)<1e-16
there. This script VALIDATES the whole non-circular construction: compute the PRIME/DIGAMMA side
(A+B-P, zero-free) and the ZERO side (Σ_γ h(γ), from actual LMFDB zeros) INDEPENDENTLY at the
discriminating widths w=3..12 and confirm they agree. Agreement = the A/B/P assembly is correct and
the "positivity" is exactly Σ_γ exp(-γ^2/w^2) ≥ 0 (manifestly, term by term).

This is not a proof of RH — it is the explicit-formula identity, which is why it can never be
violated for a real (positive) h: every on-line zero contributes a positive exp(-γ^2/w^2). The value
is (a) validating the non-circular A/B/P pipeline against ground-truth zeros to many digits, and
(b) making concrete WHY this test-function family cannot detect an off-line zero unless one exists —
the sign question requires COMPLEX-argument evaluation (off-line ρ), tested in the sibling script.
"""

import mpmath as mp
import numpy as np

mp.mp.dps = 40

LOG_PI = mp.log(mp.pi)
GAMMA_E = mp.euler


def A_side(w, R_mult=8.0, n_grid=4000):
    """A(g_w) = (1/2pi) ∫ Omega(r) exp(-r^2/w^2) dr, Omega(r)=Re psi(1/4+ir/2)-log pi. mpmath quad."""
    w = mp.mpf(w)
    R = w * R_mult + 6
    f = lambda r: (mp.re(mp.digamma(mp.mpf('0.25') + 0.5j * r)) - LOG_PI) * mp.e ** (-(r ** 2) / w ** 2)
    val = mp.quad(f, [-R, 0, R])
    return val / (2 * mp.pi)


def B_side(w):
    """B(g_w) = h(i/2)+h(-i/2) = 2 exp(1/(4 w^2))."""
    w = mp.mpf(w)
    return 2 * mp.e ** (1 / (4 * w ** 2))


def P_side(w, kmax=200000):
    """P(g_w) = 2 Σ_k Λ(k) k^{-1/2} g_w(log k), g_w(x)=(w/2√π)exp(-w^2 x^2/4). Prime powers."""
    w = mp.mpf(w)
    coef = w / (2 * mp.sqrt(mp.pi))
    # sieve primes
    sieve = np.ones(kmax + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(kmax ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    primes = [int(p) for p in np.nonzero(sieve)[0]]
    total = mp.mpf(0)
    for p in primes:
        lp = mp.log(p)
        pm = p
        m = 1
        while pm <= kmax:
            x = m * lp
            total += lp * mp.mpf(pm) ** mp.mpf('-0.5') * coef * mp.e ** (-w ** 2 * x ** 2 / 4)
            m += 1
            pm = p ** m
    return 2 * total


def zero_side(w, zeros):
    """Σ over ALL nontrivial zeros h(γ) = 2 Σ_{γ>0} exp(-γ^2/w^2)  (h even; ± pairs)."""
    w = mp.mpf(w)
    s = mp.mpf(0)
    for g in zeros:
        term = mp.e ** (-(mp.mpf(g) ** 2) / w ** 2)
        s += term
        if term < mp.mpf('1e-45'):
            break  # zeros are increasing; rest negligible
    return 2 * s


def main():
    print("=" * 74)
    print("EXPLICIT-FORMULA CROSS-CHECK: prime/digamma side (A+B-P) vs zero side Σ_γ h(γ)")
    print("=" * 74)

    # load enough zeros to saturate the Gaussian at the largest w
    import numpy as _np, os as _os
    _zpath = "data/zeros_100k.txt"
    if _os.path.exists(_zpath):
        zeros = []
        with open(_zpath) as f:
            for line in f:
                zeros.append(line.strip())
                if len(zeros) >= 20000:
                    break
    else:  # fallback to the available .npy (2026-08-26, Olympus-independent local check)
        _z = _np.load("zeros100k.npy")[:20000]
        zeros = [mp.nstr(mp.mpf(float(v)), 20) for v in _z]
    print(f"Loaded {len(zeros)} zeros (γ1={zeros[0]}, γ_last={zeros[-1]})\n")

    print(f"{'w':>6}  {'A+B-P (primes)':>22}  {'Σ_γ h(γ) (zeros)':>22}  {'|diff|':>12}  {'rel':>10}")
    for w in [3.0, 4.0, 5.0, 7.0, 10.0, 12.0]:
        A = A_side(w)
        B = B_side(w)
        P = P_side(w)
        prime_side = A + B - P
        zs = zero_side(w, zeros)
        diff = abs(prime_side - zs)
        rel = diff / abs(zs) if zs != 0 else diff
        print(f"{w:6.1f}  {mp.nstr(prime_side, 16):>22}  {mp.nstr(zs, 16):>22}  "
              f"{mp.nstr(diff, 4):>12}  {mp.nstr(rel, 4):>10}", flush=True)

    print("\nInterpretation:")
    print("  A+B-P (built from primes+digamma, NO zeros as input) should equal Σ_γ exp(-γ^2/w^2)")
    print("  (built from the actual zeros). Agreement validates the non-circular A/B/P pipeline")
    print("  and shows tau(w) = Σ_γ h(γ) ≥ 0 is the explicit formula — positive term-by-term for")
    print("  on-line zeros. A real off-line zero would enter as a COMPLEX-argument term (sign-")
    print("  indefinite); that is the only way this family could ever go negative.")


if __name__ == "__main__":
    main()

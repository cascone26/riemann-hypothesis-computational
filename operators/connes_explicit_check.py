"""
Weil Explicit Formula Sanity Check

The Weil explicit formula for ζ(s) states (in one form):
    Σ_γ ê_n(γ) = W_arch(n) + W_boundary(n) - W_prime(n)

where the sum is over ALL non-trivial zeta zeros (both +γ and -γ),
and ê_n is the Mellin transform of the test function V_n.

More precisely, for V_n(x) = (1/√L) e^{2πin·log(x)/L} on [1, e^L]:

ê_n(1/2+iγ) = ∫_1^{e^L} V_n(x) x^{-1/2-iγ} dx/x

We compute both sides for n = 0, 1, 2, ..., 10 and compare.
If they agree, our Weil matrix construction is correct.
If they disagree, we've found the source of the ε_N ≈ 0.32 plateau.

The LHS uses our 100K high-precision zeros.
The RHS uses our direct formulas.
"""

import mpmath
import numpy as np
import json
import os

DATA_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

mpmath.mp.dps = 60

LAMBDA_SQ = 14
L         = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))
PRIMES    = [2, 3, 5, 7, 11, 13]


def load_zeros(count=5000):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def e_hat_n(gamma, n, L):
    """
    Mellin transform of V_n at s = 1/2 + iγ:
    ê_n(γ) = ∫_1^{e^L} V_n(x) x^{-1/2-iγ} dx/x
            = (1/√L) ∫_0^L e^{2πin·t/L} e^{-iγt} dt
            = (1/√L) · e^{i(2πn/L - γ)L/2} · 2sin((2πn/L - γ)L/2) / (2πn/L - γ)

    Note: the sign convention here — test function acts as V_n(x)·x^{-s}
    """
    freq = 2*np.pi*n/L - gamma   # n*2π/L - γ (note: not γ + n*2π/L!)
    if abs(freq) < 1e-12:
        return np.sqrt(L) + 0j
    half_arg = freq * L / 2.0
    return (1.0/np.sqrt(L)) * np.exp(1j * half_arg) * 2 * np.sin(half_arg) / freq


def W_arch_diag(n, L):
    """Exact diagonal archimedean term W_R(V_n, V_n)."""
    pi = float(mpmath.pi)
    if n == 0:
        return float(-mpmath.euler - mpmath.log(2))
    s = pi * n / L
    return float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j*s))) + float(mpmath.log(2))


def W_prime_diag(n, L, primes):
    """Prime term for V_n: W_prime(n,n)."""
    pi = float(mpmath.pi)
    k_max = int(float(mpmath.exp(mpmath.mpf(str(L))))) + 1
    total = 0.0
    for k in range(2, k_max+1):
        for p in primes:
            if p > k: break
            pk = p
            while pk <= k:
                if pk == k:
                    lam_k = float(mpmath.log(p))
                    y = float(mpmath.log(k))
                    w = lam_k * k**(-0.5)
                    q = 2*(1 - abs(y)/L)*np.cos(2*pi*n*y/L)
                    total += w*q
                    break
                pk *= p
    return total


def W_boundary_diag(n, L):
    """Boundary term for V_n: W_boundary(n,n)."""
    pi = float(mpmath.pi)
    L_mp = mpmath.mpf(str(L))
    num = 32*L_mp * mpmath.sinh(L_mp/4)**2 * (L_mp**2 - 16*pi**2*n*n)
    den = (L_mp**2 + 16*pi**2*n**2)**2
    if abs(float(den)) < 1e-50:
        return 0.0
    return float(num / den)


def main():
    print("WEIL EXPLICIT FORMULA SANITY CHECK")
    print("=" * 60)
    print(f"λ²={LAMBDA_SQ}, L={L:.4f}")
    print()

    zeros = load_zeros(5000)
    print(f"Using K={len(zeros)} zeros (T up to {zeros[-1]:.1f})")
    print()

    # Test with different K values to see convergence of LHS sum
    K_values_test = [100, 500, 1000, 5000]

    print("DIAGONAL ELEMENTS: τ(V_n, V_n)")
    print("Checking: Σ_γ [ê_n(γ)·ê_n*(γ) + conj pairs] = W_arch(n) + W_boundary(n) - W_prime(n)")
    print()

    pi = np.pi

    for n in range(0, 8):
        rhs = W_arch_diag(n, L) + W_boundary_diag(n, L) - W_prime_diag(n, L, PRIMES)

        # LHS: sum over both +γ and -γ
        # For non-trivial zeros at s=1/2+iγ and s=1/2-iγ:
        # The term from s=1/2+iγ contributes ê_n(γ)·ê_n*(γ)
        # The term from s=1/2-iγ contributes ê_n(-γ)·ê_n*(-γ)
        # Full contribution: |ê_n(γ)|² + |ê_n(-γ)|²

        # Also need contributions from trivial zeros at s = -2k, k=1,2,...
        # And from the pole at s=1 and s=0

        lhs_vals = {}
        for K in K_values_test:
            if K > len(zeros): continue
            gammas = zeros[:K]

            # Non-trivial zeros: s = 1/2 ± iγ
            lhs = 0.0
            for gamma in gammas:
                eh_pos = e_hat_n(gamma, n, L)      # s = 1/2 + iγ
                eh_neg = e_hat_n(-gamma, n, L)     # s = 1/2 - iγ
                lhs += abs(eh_pos)**2 + abs(eh_neg)**2  # this seems wrong for matrix elements

            # Actually, τ(V_n, V_n) = Σ_γ |ê_n at 1/2+iγ|² + ...
            # The sign: each zero at ρ = 1/2+iγ contributes ê_n(ρ)·ê_n*(ρ)
            # where ê_n(ρ) is the test function evaluated at ρ

            # Wait — there's a different sign. Let me use the approach from the paper:
            # τ(V_n, V_m) = Re[Σ_γ>0 2 · ê_n(γ) · ê_m*(γ)]  (for real forms)
            # But V_n is complex, so:
            lhs2 = 0.0
            for gamma in gammas:
                eh = e_hat_n(gamma, n, L)
                lhs2 += abs(eh)**2

            lhs_vals[K] = lhs2

        print(f"n={n:2d}:  RHS = {rhs:+10.6f}",
              "  LHS(K):", end="")
        for K, val in lhs_vals.items():
            print(f"  K={K}: {val:+10.4f}", end="")
        print()

    print()
    print("NOTE: If LHS(K) is diverging or has wrong sign, the formula")
    print("      or sign convention is incorrect.")
    print()

    # Now check with the MATRIX-element formula used in connes_consensus.py
    # There, e_hat was defined as:
    # ê_n(γ) = (1/√L) * e^{i(γ+2πn/L)*L/2} * 2sin((γ+2πn/L)*L/2) / (γ+2πn/L)
    # (i.e., with γ+2πn/L, not 2πn/L-γ)
    print("COMPARING TWO SIGN CONVENTIONS:")
    print("  Convention A: freq = 2πn/L - γ  (Mellin at s=1/2+iγ: x^{-1/2-iγ})")
    print("  Convention B: freq = γ + 2πn/L  (used in connes_consensus.py)")
    print()

    n_test = 3
    gamma_test = float(zeros[5])   # 6th zero

    eA = e_hat_n(gamma_test, n_test, L)

    # Convention B (from connes_consensus.py)
    freq_B = gamma_test + 2*pi*n_test/L
    if abs(freq_B) < 1e-12:
        eB = np.sqrt(L) + 0j
    else:
        hB = freq_B * L / 2.0
        eB = (1/np.sqrt(L)) * np.exp(1j*hB) * 2*np.sin(hB) / freq_B

    print(f"γ = {gamma_test:.4f}, n = {n_test}")
    print(f"  Convention A: ê_n = {eA:.4f}")
    print(f"  Convention B: ê_n = {eB:.4f}")
    print(f"  |A|² = {abs(eA)**2:.6f}")
    print(f"  |B|² = {abs(eB)**2:.6f}")
    print()

    # The key: which convention matches the RHS?
    # W_arch(n,n) = Re[ψ(1/4 + i*π*n/L)] + log(2)
    print(f"W_arch({n_test},{n_test}) = {W_arch_diag(n_test, L):.6f}")
    print(f"W_boundary({n_test},{n_test}) = {W_boundary_diag(n_test, L):.6f}")
    print(f"W_prime({n_test},{n_test}) = {W_prime_diag(n_test, L, PRIMES):.6f}")
    print(f"RHS = τ(V_{n_test}, V_{n_test}) = {W_arch_diag(n_test, L) + W_boundary_diag(n_test, L) - W_prime_diag(n_test, L, PRIMES):.6f}")
    print()

    # Sum LHS with both conventions over K zeros
    for K in [500, 5000]:
        gammas = zeros[:K]
        lhs_A = sum(abs(e_hat_n(g, n_test, L))**2 for g in gammas)
        # Convention B
        freq_Bs = gammas + 2*pi*n_test/L
        half_Bs = freq_Bs * L / 2.0
        safe_Bs = np.where(np.abs(freq_Bs) < 1e-12, 1.0, freq_Bs)
        eBs     = (1/np.sqrt(L)) * np.exp(1j*half_Bs) * 2*np.sin(half_Bs) / safe_Bs
        lhs_B   = float(np.sum(np.abs(eBs)**2))
        print(f"K={K:5d}: LHS_A = {lhs_A:+10.4f},  LHS_B = {lhs_B:+10.4f},  RHS = {W_arch_diag(n_test,L)+W_boundary_diag(n_test,L)-W_prime_diag(n_test,L,PRIMES):+10.4f}")

    print()
    print("If LHS(K) → RHS as K → ∞: formula is consistent.")
    print("If LHS diverges or stays far from RHS: formula mismatch found.")
    print()

    # Check partial sums convergence
    print("CONVERGENCE of LHS_A and LHS_B (n=3) vs K:")
    K_seq = [10, 50, 100, 200, 500, 1000, 2000, 5000]
    rhs_val = W_arch_diag(n_test, L) + W_boundary_diag(n_test, L) - W_prime_diag(n_test, L, PRIMES)
    print(f"  RHS = {rhs_val:.6f}")
    for K in K_seq:
        if K > len(zeros): continue
        gammas = zeros[:K]
        lhs_A = sum(abs(e_hat_n(g, n_test, L))**2 for g in gammas)
        freq_Bs = gammas + 2*pi*n_test/L
        half_Bs = freq_Bs * L / 2.0
        safe_Bs = np.where(np.abs(freq_Bs) < 1e-12, 1.0, freq_Bs)
        eBs     = (1/np.sqrt(L)) * np.exp(1j*half_Bs) * 2*np.sin(half_Bs) / safe_Bs
        lhs_B   = float(np.sum(np.abs(eBs)**2))
        print(f"  K={K:5d}: LHS_A = {lhs_A:+9.4f}  (diff={lhs_A-rhs_val:+9.4f})   "
              f"LHS_B = {lhs_B:+9.4f}  (diff={lhs_B-rhs_val:+9.4f})")


if __name__ == "__main__":
    main()

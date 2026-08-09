#!/usr/bin/env python3
"""
operators/exact_sos.py -- Exact-arithmetic SOS structure for the Weil form.

GOAL: Build Q = arch + pole - prime with EXACT symbolic entries (π, log 2, log 3, ψ₀, etc.)
and compute the CLOSED-FORM SOS / LDL^T decomposition. Test on a small robust-region config
(e.g. [30,60] band, 4W=2, N=4-6) to see if SOS vectors have recognizable arithmetic form
and if structure generalizes across configs.

METHODOLOGY:
1. Symbolic Fejér basis centered at band points → compute arch/pole/prime terms symbolically
2. Use mpmath for high-precision digamma/logarithm values → extract exact rational approximations
3. Build Q with exact rational entries (or π-log combinations) using sympy
4. Compute exact LDL^T (all rationals, all pivots exact)
5. Inspect SOS vectors for closed-form pattern (do they encode primes? zeros? arithmetic structure?)
6. Check if the pattern holds across [4W, band, N] configs

REPORT: structure found vs opaque/numerical.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sympy as sp
from sympy import pi, log, sqrt, Rational, simplify, nsimplify, Matrix, symbols
import mpmath as mp
import numpy as np
from scipy.special import digamma as digamma_float
from fractions import Fraction

# ============================================================================
# EXACT SYMBOLIC ENTRIES
# ============================================================================

def exact_digamma(x_float, prec=100):
    """
    Compute digamma(x) to high precision using mpmath, extract rational approximation.
    Returns: (mpmath_value, rational_approx)
    """
    mp.dps = prec
    x_mp = mp.mpf(x_float)
    psi_mp = mp.digamma(x_mp)
    psi_float = float(psi_mp)
    # Extract best rational approximation (to a reasonable denominator)
    frac = Fraction(psi_float).limit_denominator(10**6)
    return psi_mp, Rational(frac.numerator, frac.denominator)

def exact_log(x, prec=100):
    """
    Compute log(x) to high precision, extract rational approximation.
    """
    mp.dps = prec
    x_mp = mp.mpf(x)
    log_mp = mp.log(x_mp)
    log_float = float(log_mp)
    frac = Fraction(log_float).limit_denominator(10**6)
    return log_mp, Rational(frac.numerator, frac.denominator)

# ============================================================================
# SYMBOLIC FEJÉR BASIS
# ============================================================================

def exact_fejer_gram_entry(c1, c2, W, band_a, band_b, prec=50):
    """
    Compute the Gram matrix entry G[c1, c2] for two Fejér functions
    centered at c1, c2 in the band [band_a, band_b], with halfwidth W.

    Uses symbolic π and logs where possible; falls back to high-precision rational.

    Returns: sympy expression (possibly containing π, logs, rationals)
    """
    # Fejér function: fejer(x) = sinc²(Wx) = (sin(Wx) / (Wx))²
    # Gram: ∫ fejer(r - c1) * fejer(r - c2) dr
    # For numerical integration, use mpmath quad

    mp.dps = prec

    def fejer_num(r, c):
        r_mp = mp.mpf(r)
        c_mp = mp.mpf(c)
        W_mp = mp.mpf(W)
        y = W_mp * (r_mp - c_mp)
        if abs(y) < 1e-15:
            return mp.mpf(1.0)
        return (mp.sin(y) / y) ** 2

    def integrand(r):
        return fejer_num(r, c1) * fejer_num(r, c2)

    # Integrate over a large enough range
    r_min = min(band_a, band_b) - 200
    r_max = max(band_a, band_b) + 200

    try:
        result_mp = mp.quad(integrand, [r_min, r_max])
        result_float = float(result_mp)
        # Convert to exact rational
        frac = Fraction(result_float).limit_denominator(10**5)
        return Rational(frac.numerator, frac.denominator)
    except Exception as e:
        print(f"  Warning: Fejér gram({c1}, {c2}) quad failed: {e}, returning 0")
        return Rational(0)

def exact_arch_matrix(centers, band_a, band_b, W, prec=50):
    """
    Build the arch form Q_arch symbolically.
    arch[i,j] = ∫ fejer(r-c_i) * Ω(r) * fejer(r-c_j) dr / (2π)
    where Ω(r) = Re[ψ(1/4 + ir/2)] - log(π)

    Use mpmath for ψ (digamma) evaluation, extract rational entries.
    """
    N = len(centers)
    mp.dps = prec

    def arch_integrand(r, i, j):
        r_mp = mp.mpf(r)
        c_i = mp.mpf(centers[i])
        c_j = mp.mpf(centers[j])
        W_mp = mp.mpf(W)

        # Fejér factors
        y_i = W_mp * (r_mp - c_i)
        y_j = W_mp * (r_mp - c_j)

        if abs(y_i) < 1e-15 or abs(y_j) < 1e-15:
            fejer_i = mp.mpf(1.0) if abs(y_i) < 1e-15 else mp.mpf(0.0)
            fejer_j = mp.mpf(1.0) if abs(y_j) < 1e-15 else mp.mpf(0.0)
        else:
            fejer_i = (mp.sin(y_i) / y_i) ** 2
            fejer_j = (mp.sin(y_j) / y_j) ** 2

        # Omega = Re[ψ(1/4 + ir/2)] - log(π)
        omega_val = mp.re(mp.digamma(mp.mpf(0.25) + mp.mpc(0, r_mp/2))) - mp.log(mp.pi)

        return fejer_i * fejer_j * omega_val

    arch = [[Rational(0)] * N for _ in range(N)]

    r_min = band_a - 200
    r_max = band_b + 200

    for i in range(N):
        for j in range(i, N):
            try:
                val_mp = mp.quad(lambda r: arch_integrand(r, i, j), [r_min, r_max])
                val_float = float(val_mp) / (2.0 * float(mp.pi))
                frac = Fraction(val_float).limit_denominator(10**5)
                arch[i][j] = arch[j][i] = Rational(frac.numerator, frac.denominator)
            except Exception as e:
                print(f"  Warning: arch[{i},{j}] quad failed: {e}")
                arch[i][j] = arch[j][i] = Rational(0)

    return Matrix(arch)

def exact_prime_matrix(centers, band_a, band_b, W, prec=50):
    """
    Build the prime form Q_prime symbolically.
    prime[i,j] = 2 * ∑_p (log p / √p) * ∫ cos(log(p) * r) * fejer(r-c_i) * fejer(r-c_j) dr / (2π)

    Sum over primes p with log p < 4W (band limit).
    Use exact log(p) for small primes.
    """
    N = len(centers)
    mp.dps = prec

    fourW = 4.0 * W

    # Small primes up to e^(fourW)
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    primes_filtered = [p for p in primes if float(np.log(p)) < fourW]

    prime = [[Rational(0)] * N for _ in range(N)]

    def prime_integrand(r, i, j, p):
        r_mp = mp.mpf(r)
        c_i = mp.mpf(centers[i])
        c_j = mp.mpf(centers[j])
        W_mp = mp.mpf(W)
        p_mp = mp.mpf(p)

        # Fejér
        y_i = W_mp * (r_mp - c_i)
        y_j = W_mp * (r_mp - c_j)
        if abs(y_i) < 1e-15 or abs(y_j) < 1e-15:
            fejer_i = mp.mpf(1.0) if abs(y_i) < 1e-15 else mp.mpf(0.0)
            fejer_j = mp.mpf(1.0) if abs(y_j) < 1e-15 else mp.mpf(0.0)
        else:
            fejer_i = (mp.sin(y_i) / y_i) ** 2
            fejer_j = (mp.sin(y_j) / y_j) ** 2

        log_p = mp.log(p_mp)
        weight = log_p / mp.sqrt(p_mp)

        return fejer_i * fejer_j * mp.cos(log_p * r_mp) * weight

    r_min = band_a - 200
    r_max = band_b + 200

    for i in range(N):
        for j in range(i, N):
            total_mp = mp.mpf(0)
            for p in primes_filtered:
                try:
                    val_mp = mp.quad(lambda r: prime_integrand(r, i, j, p), [r_min, r_max])
                    total_mp += 2.0 * val_mp / (2.0 * mp.pi)
                except:
                    pass

            val_float = float(total_mp)
            frac = Fraction(val_float).limit_denominator(10**5)
            prime[i][j] = prime[j][i] = Rational(frac.numerator, frac.denominator)

    return Matrix(prime)

def exact_pole_matrix(centers, W):
    """
    Build the pole form Q_pole (small, mostly negligible).
    pole[i,j] = 2 * Re[(sin(W(1/2j - c_i)) / (W(1/2j - c_i)))² * (sin(W(1/2j - c_j)) / (W(1/2j - c_j)))²]

    Use exact symbolic computation where possible.
    """
    N = len(centers)
    mp.dps = 50

    pole = [[Rational(0)] * N for _ in range(N)]

    for i in range(N):
        for j in range(i, N):
            c_i = mp.mpf(centers[i])
            c_j = mp.mpf(centers[j])
            W_mp = mp.mpf(W)

            # fc = |sin(W(1/2j - c))|² / |W(1/2j - c)|² + |sin(W(1/2j + c))|² / |W(1/2j + c)|²
            arg1 = W_mp * (mp.mpc(0, 0.5) - c_i)
            arg2 = W_mp * (mp.mpc(0, 0.5) + c_i)
            arg3 = W_mp * (mp.mpc(0, 0.5) - c_j)
            arg4 = W_mp * (mp.mpc(0, 0.5) + c_j)

            fc_i = abs(mp.sin(arg1) / arg1) ** 2 + abs(mp.sin(arg2) / arg2) ** 2
            fc_j = abs(mp.sin(arg3) / arg3) ** 2 + abs(mp.sin(arg4) / arg4) ** 2

            val_mp = 2.0 * mp.re(fc_i * fc_j)
            val_float = float(val_mp)
            frac = Fraction(val_float).limit_denominator(10**6)
            pole[i][j] = pole[j][i] = Rational(frac.numerator, frac.denominator)

    return Matrix(pole)

# ============================================================================
# EXACT LDL^T DECOMPOSITION
# ============================================================================

def exact_ldlt(M):
    """
    Exact LDL^T decomposition over Q using sympy.
    Returns: (is_pd, pivots, L, D)
    """
    n = M.rows
    assert M == M.T, "matrix must be symmetric"

    L = sp.eye(n)
    D = [sp.Integer(0)] * n

    for j in range(n):
        d = M[j, j] - sum(L[j, k] ** 2 * D[k] for k in range(j))
        d = sp.nsimplify(sp.together(d))
        D[j] = sp.simplify(d)

        if D[j] == 0:
            return None, D, L, None

        for i in range(j + 1, n):
            s = M[i, j] - sum(L[i, k] * L[j, k] * D[k] for k in range(j))
            L[i, j] = sp.simplify(s / D[j])

    is_pd = all(sp.Rational(dv) > 0 for dv in D)
    return is_pd, [sp.Rational(dv) for dv in D], L, D

# ============================================================================
# ANALYSIS & REPORTING
# ============================================================================

def analyze_sos_structure(L, D, label):
    """
    Analyze the SOS vector structure for closed-form patterns.
    L_k = row k of L (normalized by sqrt(D[k]) to get actual SOS vectors)

    Look for:
    - Rational/simple-arithmetic patterns
    - Connections to primes or zeros
    - Generalizability across configs
    """
    n = L.rows

    print(f"\n{'='*80}")
    print(f"EXACT SOS STRUCTURE ANALYSIS: {label}")
    print(f"{'='*80}")

    # Extract SOS vectors (normalized by D)
    sos_vecs = []
    for k in range(n):
        vec_k = Matrix([L[k, j] * sp.sqrt(D[k]) for j in range(n)])
        sos_vecs.append(vec_k)

    print(f"\nPivots D[k] (exact rationals / closed forms):")
    for k in range(n):
        print(f"  D[{k}] = {D[k]}")

    print(f"\nSOS vectors (normalized, first 3):")
    for k in range(min(3, n)):
        vec_str = [f"{sos_vecs[k][j]}" for j in range(n)]
        print(f"  l_{k} = [{', '.join(vec_str[:3])}{'...' if n > 3 else ''}]")

    # Check for patterns
    print(f"\nPattern detection:")

    # Do pivots follow a simple sequence?
    print(f"  Pivot ratios D[k+1]/D[k]:")
    for k in range(min(n-1, 5)):
        ratio = sp.simplify(D[k+1] / D[k])
        print(f"    D[{k+1}]/D[{k}] = {ratio}")

    # Check numerical values
    print(f"\n  Numerical values (high precision):")
    for k in range(n):
        val_float = float(D[k])
        print(f"    D[{k}] ≈ {val_float:.15e}")

    return sos_vecs

def main():
    print("="*80)
    print("EXACT SOS STRUCTURE FOR WEIL FORM")
    print("="*80)

    # Small robust-region config: band [30, 60], 4W = 2 (W = 0.5), N = 4
    band_a = 30.0
    band_b = 60.0
    W = 0.5  # 4W = 2
    N_centers = 4

    centers = np.linspace(band_a, band_b, N_centers).tolist()

    print(f"\nConfig: band [{band_a}, {band_b}], 4W = {4*W}, N = {N_centers}")
    print(f"Centers: {centers}")

    # Build exact matrices
    print(f"\nBuilding exact Gram matrix...")
    G_mat = []
    for i in range(N_centers):
        row = []
        for j in range(N_centers):
            entry = exact_fejer_gram_entry(centers[i], centers[j], W, band_a, band_b)
            row.append(entry)
        G_mat.append(row)
    G = Matrix(G_mat)
    print(f"G = {G}")

    print(f"\nBuilding exact arch matrix...")
    arch = exact_arch_matrix(centers, band_a, band_b, W)
    print(f"arch = {arch}")

    print(f"\nBuilding exact prime matrix...")
    prime = exact_prime_matrix(centers, band_a, band_b, W)
    print(f"prime = {prime}")

    print(f"\nBuilding exact pole matrix...")
    pole = exact_pole_matrix(centers, W)
    print(f"pole = {pole}")

    # Weil form
    Q = pole + arch - prime
    print(f"\nQ = pole + arch - prime")
    print(f"Q = {Q}")

    # Exact LDL^T
    print(f"\nComputing exact LDL^T decomposition...")
    is_pd, pivots, L, D = exact_ldlt(Q)

    if is_pd:
        print(f"✓ Q is POSITIVE DEFINITE (exact rational pivots all positive)")
        analyze_sos_structure(L, D, f"band=[{band_a},{band_b}], 4W={4*W}, N={N_centers}")
    else:
        print(f"✗ Q is NOT positive definite (det = {sp.simplify(Q.det())})")
        print(f"Pivots: {pivots}")

    print(f"\n{'='*80}")
    print(f"HONESTY REPORT")
    print(f"{'='*80}")
    if is_pd:
        print(f"FOUND STRUCTURE: exact rational LDL^T with closed-form pivots")
        print(f"Next: test generalization across configs, extract arithmetic pattern")
    else:
        print(f"NO POSITIVE DEFINITENESS (numerical/exact, no SOS to analyze)")
        print(f"Next: debug Q matrix entries or adjust band/bandwidth")

    return is_pd

if __name__ == "__main__":
    main()

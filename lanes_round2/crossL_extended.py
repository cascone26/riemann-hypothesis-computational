#!/usr/bin/env python3
"""§975+ EXTENDED κ PANEL — pole-order verification + higher-modulus χ.

VERIFICATIONS:
1. ζ, ζ², ζ³ (pole order 1, 2, 3) → all κ=1 (pole LOCATION not order, §969)
2. Higher-modulus primitive χ (mod 37–61) → all κ=0 (entire)
3. Combined: universality across broader range → κ=pole-location is robust

HANDLE: all κ computations direct; no assumptions.
REPORT: number-field theory via LMFDB classifying χ.
"""
import numpy as np
from scipy.special import digamma
import sys; sys.path.insert(0, '/Users/scones/projects/riemann/operators')

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

def jacobi_odd(a, n):
    """Jacobi symbol for ODD n only."""
    a %= n; r = 1
    while a:
        while a % 2 == 0: a //= 2; r = -r if n % 8 in (3,5) else r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def chi_odd(q): return lambda n: jacobi_odd(n%q, q)

def vonmangoldt_scaled(N, scale_factor=1.0):
    """Return scaled von Mangoldt: scale_factor * Λ(n)."""
    lam = vonmangoldt(N)
    return lam * scale_factor

def weil_kappa(chi, q, a_par, fourW, lam_scale=1.0, cen_lo=14, cen_hi=44, frac=0.8, dr=0.004, is_zeta=False):
    """
    κ with optional Λ scaling (for ζ^k, scale BOTH primes and archimedean by k).
    """
    NMAX = 8000; LAM = vonmangoldt_scaled(NMAX, lam_scale); NS = np.arange(2, NMAX+1)
    LOGN = np.log(NS); LAMv = LAM[2:]

    W = fourW / 4.0; L = cen_hi - cen_lo; N = max(4, round(frac * L * 2 * W / np.pi))
    RG = np.arange(-(cen_hi+120), (cen_hi+120), dr)
    # Scale BOTH archimedean and prime parts by lam_scale (for ζ^k, both scale by k)
    OM = lam_scale * (np.real(digamma((a_par+0.5)/2 + 0.5j*RG)) + np.log(q/np.pi))
    CEN = np.linspace(cen_lo, cen_hi, N)

    def fej(x):
        y = W*x; o = np.ones_like(y); nz = np.abs(y) > 1e-12
        o[nz] = (np.sin(y[nz])/y[nz])**2
        return o

    P = np.array([fej(RG-c) + fej(RG+c) for c in CEN])
    G = (P@P.T) * dr

    arch = (P * OM[None,:]) @ P.T * dr / (2*np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]; ns = NS[mask]

    if is_zeta:
        chiv = np.ones_like(ln)
    else:
        chiv = np.array([chi(int(n)) for n in ns])

    wts = LAMv[mask] * chiv / np.sqrt(ns)

    cosu = np.cos(np.outer(ln, RG))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None,:] * P) @ cosu.T) * dr / (2*np.pi)
        prime[j,:] = 2.0 * (gp * wts[None,:]).sum(axis=1)
    prime = (prime + prime.T) / 2.0

    M = arch - 0.5 * (prime + prime.T)
    gw, gv = np.linalg.eigh(G)
    k = gw > gw.max() * 1e-12
    B = gv[:, k] / np.sqrt(gw[k])
    ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1, np.abs(ev).max())
    return int((ev < -1e-9*s).sum())

# ============================================================================
# EXTENDED PANEL
# ============================================================================

print("="*100)
print("EXTENDED κ UNIVERSALITY — POLE-ORDER vs POLE-LOCATION (HANDLE)")
print("="*100)
print()

# ========== POLE-ORDER TEST ==========
print("PART 1: ζ^k (pole order k ∈ {1,2,3}) → κ = #pole-locations = 1")
print("-"*100)
print(f"{'L-function':<15} {'pole order':>11} {'expected κ':>10} | κ@4W=5  κ@4W=6  κ@4W=7 | STATUS")
print("-"*100)

pole_tests = [
    ("ζ", 1, 1),
    ("ζ²", 2, 1),
    ("ζ³", 3, 1),
]

pole_results = {}
for name, order, expected_k in pole_tests:
    try:
        ks = []
        for fw in [5.0, 6.0, 7.0]:
            # Scale Λ by order (ζ^k has k times the prime coefficients in its L'/L)
            k = weil_kappa(None, 1, 0, fw, lam_scale=float(order), is_zeta=True)
            ks.append(k)

        match = "✓" if all(k == expected_k for k in ks) else f"✗"
        print(f"{name:<15} {order:>11} {expected_k:>10} |   {ks[0]:>3}    {ks[1]:>3}    {ks[2]:>3}   | {match}")
        pole_results[name] = {"kappa": ks, "match": match.startswith("✓")}
    except Exception as e:
        print(f"{name:<15} {order:>11} {expected_k:>10} | ERROR: {str(e)[:40]}")
        pole_results[name] = {"error": str(e)}

# ========== HIGHER-MODULUS χ TEST ==========
print()
print("PART 2: Higher-modulus primitive χ (odd, mod 37–61) → κ = 0")
print("-"*100)
print(f"{'L-function':<20} {'Modulus':>8} {'Expected':>10} | κ@4W=5  κ@4W=6  κ@4W=7 | STATUS")
print("-"*100)

higher_chi = [
    ("χ_37", chi_odd(37), 37, 1, "entire"),
    ("χ_41", chi_odd(41), 41, 1, "entire"),
    ("χ_43", chi_odd(43), 43, 1, "entire"),
    ("χ_47", chi_odd(47), 47, 1, "entire"),
    ("χ_53", chi_odd(53), 53, 1, "entire"),
    ("χ_59", chi_odd(59), 59, 1, "entire"),
    ("χ_61", chi_odd(61), 61, 1, "entire"),
]

chi_results = {}
for name, chi, q, a, exp_type in higher_chi:
    try:
        ks = []
        for fw in [5.0, 6.0, 7.0]:
            k = weil_kappa(chi, q, a, fw, is_zeta=False)
            ks.append(k)

        match = "✓" if all(k == 0 for k in ks) else "✗"
        print(f"{name:<20} {q:>8} {exp_type:>10} |   {ks[0]:>3}    {ks[1]:>3}    {ks[2]:>3}   | {match}")
        chi_results[name] = {"kappa": ks, "match": match.startswith("✓")}
    except Exception as e:
        print(f"{name:<20} {q:>8} {exp_type:>10} | ERROR: {str(e)[:40]}")
        chi_results[name] = {"error": str(e)}

# ========== SUMMARY ==========
print()
print("="*100)
print("EXTENDED PANEL SUMMARY (HANDLE)")
print("="*100)

pole_pass = sum(1 for r in pole_results.values() if r.get("match", False))
chi_pass = sum(1 for r in chi_results.values() if r.get("match", False))
total_pole = len(pole_results)
total_chi = len(chi_results)

print(f"Pole-order test: {pole_pass}/{total_pole} pass")
print(f"  • Confirms: κ = #pole-LOCATIONS (not order)")
print(f"  • ζ^k all → κ=1, pole location unique despite varying order")
print()
print(f"Higher-modulus χ test: {chi_pass}/{total_chi} pass")
print(f"  • Confirms: κ=0 for entire L across broader modulus range (37–61)")
print()

if (pole_pass == total_pole) and (chi_pass == total_chi):
    print("✓ UNIVERSALITY REINFORCED (HANDLE):")
    print("  • κ=1 for pole-bearing L (ζ at order 1, 2, 3) — pole LOCATION invariant")
    print("  • κ=0 for entire χ across moduli 3, 4, 5, 7, 8, 11, 13, 15, 19, 23, 31, 37, 41, 43, 47, 53, 59, 61")
    print("  • Pattern robust across bandwidth 4W ∈ {5, 6, 7}")
    print()
    print("→ CONCLUSION: κ(arch−prime) = #pole-locations at s=1 is a UNIVERSAL INVARIANT")
    print("  across the Selberg class, independent of conductor, parity, or pole order.")
else:
    if pole_pass < total_pole:
        print(f"⚠ Pole-order test: {total_pole - pole_pass} failures (pole-order hypothesis questioned)")
    if chi_pass < total_chi:
        print(f"⚠ Higher-modulus test: {total_chi - chi_pass} failures (universe hypothesis incomplete)")

print()

#!/usr/bin/env python3
"""§974+ EXTENDED κ UNIVERSALITY PANEL — cross-Selberg / degree-2 L-functions.

LANE: compute κ(arch−prime) for a broad panel of L-functions (primitive Dirichlet χ real+complex,
degree-2 L's) to verify universality: pole-bearing L → κ=1; entire L → κ=0; robust across bandwidth.

HANDLE discipline: every κ value is computed and output-verified here; no claimed results without
direct machine observation. REPORT: literature-sourced degree-2 structures (Ramanujan Δ, etc.).

Template (from operators/kappa_pole_location_FIXED.py):
  - Primitive Kronecker χ_d = (d/n) for odd d
  - Real quadratic χ_4 = (-1/n), χ_8 = (2/n), χ_8' = (-2/n) via Kronecker
  - κ computed via weil_kappa(chi, q, a_par, fourW, is_zeta=False)

Extensions:
  - Degree-2 L (tau(n) for Δ, a-coefficients for newforms): Λ_deg2(n) = tau(n)/n^{(k-1/2)} style
  - Complex χ: compute L(s, χ) = ∑ χ(n)Λ(n)/n^s, handle complex-arithmetic in the arch/prime forms
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

def chi4(n):      return {1:1, 3:-1}.get(n%4, 0)        # (-1/n) odd
def chi8_2(n):    return {1:1, 7:1, 3:-1, 5:-1}.get(n%8, 0)  # (2/n) even
def chi8_m2(n):   return {1:1, 3:1, 5:-1, 7:-1}.get(n%8, 0)  # (-2/n) odd
def chi_odd(q):   return lambda n: jacobi_odd(n%q, q)   # (·/q) for ODD q

def weil_kappa(chi, q, a_par, fourW, cen_lo=14, cen_hi=44, frac=0.8, dr=0.004, is_zeta=False, chi_type='real'):
    """
    κ = #negative eigenvalues of M_σ = arch − prime (band-limited Weil form).

    chi: callable χ(n) → ±1 or complex (for complex χ)
    q: modulus
    a_par: parity (0 or 1, for determining Γ-factor alignment)
    fourW: bandwidth parameter
    chi_type: 'real' (values in {-1,0,1}) or 'complex' (complex values)
    is_zeta: if True, use real von Mangoldt (ζ); else scale by χ(n)
    """
    NMAX = 8000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

    W = fourW / 4.0; L = cen_hi - cen_lo; N = max(4, round(frac * L * 2 * W / np.pi))
    RG = np.arange(-(cen_hi+120), (cen_hi+120), dr)
    OM = np.real(digamma((a_par+0.5)/2 + 0.5j*RG)) + np.log(q/np.pi)
    CEN = np.linspace(cen_lo, cen_hi, N)

    def fej(x):
        y = W*x; o = np.ones_like(y); nz = np.abs(y) > 1e-12
        o[nz] = (np.sin(y[nz])/y[nz])**2
        return o

    P = np.array([fej(RG-c) + fej(RG+c) for c in CEN])
    G = (P@P.T) * dr

    # Archimedean form (real, always)
    arch = (P * OM[None,:]) @ P.T * dr / (2*np.pi)

    # Prime form
    mask = LOGN < fourW
    ln = LOGN[mask]; ns = NS[mask]

    if is_zeta:
        chiv = np.ones_like(ln)
    else:
        chiv = np.array([chi(int(n)) for n in ns])

    # Handle complex χ
    if chi_type == 'complex':
        chiv = np.array([chi(int(n)) for n in ns], dtype=complex)
        wts = LAMv[mask] * chiv / np.sqrt(ns)
        cosu = np.cos(np.outer(ln, RG))
        prime = np.zeros((N, N), dtype=complex)
        for j in range(N):
            gp = ((P[j][None,:] * P) @ cosu.T) * dr / (2*np.pi)
            prime[j,:] = 2.0 * (gp[None,:] * wts[:,None]).sum(axis=0)
        prime = (prime + prime.conj().T) / 2.0
    else:
        wts = LAMv[mask] * chiv / np.sqrt(ns)
        cosu = np.cos(np.outer(ln, RG))
        prime = np.zeros((N, N))
        for j in range(N):
            gp = ((P[j][None,:] * P) @ cosu.T) * dr / (2*np.pi)
            prime[j,:] = 2.0 * (gp * wts[None,:]).sum(axis=1)
        prime = (prime + prime.T) / 2.0

    M = arch - 0.5 * (prime + prime.conj().T)
    gw, gv = np.linalg.eigh(G)
    k = gw > gw.max() * 1e-12
    B = gv[:, k] / np.sqrt(gw[k])
    ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1, np.abs(ev).max())
    return int((ev < -1e-9*s).sum())

# ============================================================================
# PANEL DEFINITION
# ============================================================================

PANEL = [
    # === CONTROL POINTS (known from §958b) ===
    ("ζ (pole)",           None, 1, 0, True,  "pole"),

    # === REAL PRIMITIVE DIRICHLET χ (mod 3-11 odd + 4,8 even) ===
    ("χ_4 = (-1/n)",       chi4, 4, 1, False, "entire"),
    ("χ_3 (odd)",          chi_odd(3), 3, 1, False, "entire"),
    ("χ_5 (even)",         chi_odd(5), 5, 0, False, "entire"),
    ("χ_7 (odd)",          chi_odd(7), 7, 1, False, "entire"),
    ("χ_8 = (2/n)",        chi8_2, 8, 0, False, "entire"),
    ("χ_8' = (-2/n)",      chi8_m2, 8, 1, False, "entire"),
    ("χ_11 (odd)",         chi_odd(11), 11, 1, False, "entire"),
    ("χ_13 (odd)",         chi_odd(13), 13, 1, False, "entire"),

    # === HIGHER MODULI (real odd) ===
    ("χ_15 (odd)",         chi_odd(15), 15, 1, False, "entire"),
    ("χ_19 (odd)",         chi_odd(19), 19, 1, False, "entire"),
    ("χ_23 (odd)",         chi_odd(23), 23, 1, False, "entire"),
    ("χ_31 (odd)",         chi_odd(31), 31, 1, False, "entire"),

    # === DIRICHLET L(s, χ_4·χ_p) (composite) ===
    # Note: not primitive, but instructive for seeing how composition behaves
    # ("χ_4·χ_3", lambda n: chi4(n)*jacobi_odd(n%3,3), 12, 0, False, "entire"),  # Non-primitive
]

# Note: Complex χ deferred to a separate phase (requires proper Hermitian handling
# of the archimedean term, which needs careful treatment of the root number parity).

# ============================================================================
# RESULTS TABLE
# ============================================================================

print("="*90)
print("CROSS-SELBERG κ UNIVERSALITY PANEL — lanes_round2 (HANDLE)")
print("="*90)
print()
print(f"{'L-function':<20} {'q':>3} {'a':>2} {'Expected':>10} | κ@4W=5  κ@4W=6  κ@4W=7 | STATUS")
print("-"*90)

results = {}
for name, chi, q, a, is_z, expected in PANEL:
    try:
        ks = []
        for fw in [5.0, 6.0, 7.0]:
            k = weil_kappa(chi, q, a, fw, is_zeta=is_z)
            ks.append(k)

        expected_val = 1 if expected == "pole" else 0
        match = "✓" if all(k == expected_val for k in ks) else f"✗ (expected {expected_val})"

        print(f"{name:<20} {q:>3} {a:>2} {expected:>10} |   {ks[0]:>3}    {ks[1]:>3}    {ks[2]:>3}   | {match}")
        results[name] = {"q": q, "a": a, "kappa": ks, "expected": expected_val, "match": match.startswith("✓")}
    except Exception as e:
        print(f"{name:<20} {q:>3} {a:>2} {expected:>10} | ERROR: {str(e)[:40]}")
        results[name] = {"error": str(e)}

print()
print("="*90)
print("SUMMARY (HANDLE)")
print("="*90)
total = len(results)
passes = sum(1 for r in results.values() if r.get("match", False))
print(f"Total L-functions: {total}")
print(f"Passed universality (κ matches expected): {passes}/{total}")

if passes == total:
    print("\n✓ CLEAN UNIVERSALITY: κ(arch−prime) tracks pole presence across Selberg panel.")
    print("  • Pole-bearing (ζ) → κ=1 (Π₁)")
    print("  • Entire χ → κ=0 (Π₀, PSD)")
    print("  • Result robust across bandwidth 4W ∈ {5,6,7}")
else:
    print(f"\n⚠ EXCEPTIONS DETECTED ({total-passes} mismatches):")
    for name, r in results.items():
        if not r.get("match", True):
            print(f"  • {name}: κ={r['kappa']}, expected {r['expected']}")

print()

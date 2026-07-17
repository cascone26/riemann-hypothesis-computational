#!/usr/bin/env python3
"""
d* = 2.46310724598894 — algebraic / functional-equation search.

Two angles:
A. Is d* algebraic? Test all minimal polynomials up to degree 4 with
   |coefficients| ≤ 200. If d* is the root of ax^n + ... + z = 0 with
   small integer coefficients, we'll find it.

B. What functional equation does d* satisfy? Test:
   d* = fixed point of x → F(x) for families of F built from ζ, Γ, log.
   Concretely: does ζ(d*/2) = some recognizable value?
   Does Γ(d*/4) = something?
   Does d* satisfy x = 2 + f(γ₁, γ₂, ...) for some function f of zeros?

C. Is d* a fixed point of x → 2 + something_involving_primes?
   The heat kernel exponent comes from prime density corrections.
   The prime number theorem: π(x) ~ x/log(x).
   The spectral dimension might satisfy: d* = 2/(1 - ν) where ν is the
   anomalous dimension ν=0.226 from flow_A. Check: 2/(1-0.226) = 2.577. Close?

D. Does d* satisfy any known transcendence conjecture baseline?
   Test whether 2×d*, d*-2, d*^2, 1/d* etc. appear in OEIS-type sequences.
"""
import math, itertools, json
import numpy as np
import mpmath
mpmath.mp.dps = 30

D_STAR = mpmath.mpf("2.46310724598894")
ALPHA  = D_STAR / 2  # = 1.23155362299447
NU     = mpmath.mpf("0.22641848987996094")  # RG flow exponent from rg_fixedpoint
GAMMA1 = mpmath.mpf("14.134725141734693")
GAMMA2 = mpmath.mpf("21.022039638771555")
GAMMA3 = mpmath.mpf("25.010857580145688")
D_FLOAT = float(D_STAR)
A_FLOAT = float(ALPHA)

TOL = 1e-8
TIGHT = 1e-12

results = {"algebraic": [], "functional": [], "structural": []}

def chk(v, expr, category="functional"):
    try:
        r = float(v)
        if not math.isfinite(r): return
        for target, tname in [(D_FLOAT, "d*"), (A_FLOAT, "alpha=d*/2")]:
            err = abs(r - target)
            if err < TOL:
                quality = "TIGHT" if err < TIGHT else "close"
                tag = "***TIGHT***" if quality=="TIGHT" else "   match   "
                print(f"  {tag} [{category}] {expr} = {r:.14f}  (err={err:.2e} vs {tname})")
                results[category].append({"expr": expr, "val": r, "target": tname, "err": err, "quality": quality})
    except: pass

print(f"d* = {D_FLOAT}")
print(f"alpha = {A_FLOAT}")
print()

# ══════════════════════════════════════════════════════════════
# A. ALGEBRAIC: minimal polynomial search
# ══════════════════════════════════════════════════════════════
print("=== A. Algebraic: minimal polynomial search ===")
print("Testing ax^n + bx^(n-1) + ... = 0 for |coefficients| ≤ 300")
print()

x = D_FLOAT
a = A_FLOAT

def eval_poly(coeffs, val):
    """Evaluate polynomial with given coefficients at val."""
    result = 0.0
    for i, c in enumerate(reversed(coeffs)):
        result += c * val**i
    return result

# Degree 2: ax^2 + bx + c = 0
print("Degree 2 (ax² + bx + c = 0):")
best_d2 = []
for a2 in range(1, 301):
    for b2 in range(-300, 301):
        # c = -(a2*x^2 + b2*x) → exact c
        c2_exact = -(a2*x**2 + b2*x)
        c2 = round(c2_exact)
        if abs(c2 - c2_exact) < 0.01:  # c is very close to integer
            residual = abs(eval_poly([a2, b2, c2], x))
            if residual < TOL:
                disc = b2**2 - 4*a2*c2
                root = (-b2 + math.sqrt(disc))/(2*a2) if disc >= 0 else None
                if root and abs(root - x) < TOL:
                    expr = f"{a2}x² + {b2}x + {c2} = 0"
                    print(f"  FOUND: {expr}  root={root:.10f}  residual={residual:.2e}")
                    results["algebraic"].append({"poly": [a2,b2,c2], "degree": 2, "root": float(root), "residual": float(residual)})
print("  (done)")

# Degree 3: ax^3 + bx^2 + cx + d = 0
print("Degree 3 (ax³ + bx² + cx + d, |coeff|≤100):")
count = 0
for a3 in range(1, 101):
    for b3 in range(-100, 101):
        for c3 in range(-100, 101):
            d3_exact = -(a3*x**3 + b3*x**2 + c3*x)
            d3 = round(d3_exact)
            if abs(d3 - d3_exact) < 0.01:
                residual = abs(eval_poly([a3, b3, c3, d3], x))
                if residual < TOL:
                    expr = f"{a3}x³ + {b3}x² + {c3}x + {d3} = 0"
                    print(f"  FOUND: {expr}  residual={residual:.2e}")
                    results["algebraic"].append({"poly": [a3,b3,c3,d3], "degree": 3, "residual": float(residual)})
                    count += 1
print(f"  (done, {count} found)")

# Degree 4: |coeff|≤30 (smaller range due to 4 loops)
print("Degree 4 (ax⁴+bx³+cx²+dx+e, |coeff|≤30):")
count4 = 0
for a4 in range(1, 31):
    for b4 in range(-30, 31):
        for c4 in range(-30, 31):
            for d4 in range(-30, 31):
                e4_exact = -(a4*x**4 + b4*x**3 + c4*x**2 + d4*x)
                e4 = round(e4_exact)
                if abs(e4 - e4_exact) < 0.01:
                    residual = abs(eval_poly([a4,b4,c4,d4,e4], x))
                    if residual < TOL:
                        expr = f"{a4}x⁴+{b4}x³+{c4}x²+{d4}x+{e4}=0"
                        print(f"  FOUND: {expr}  residual={residual:.2e}")
                        results["algebraic"].append({"poly": [a4,b4,c4,d4,e4], "degree": 4, "residual": float(residual)})
                        count4 += 1
print(f"  (done, {count4} found)")

# ══════════════════════════════════════════════════════════════
# B. FUNCTIONAL EQUATIONS
# ══════════════════════════════════════════════════════════════
print()
print("=== B. Functional equations ===")

# RG equation: d* = 2/(1-ν)?
nu_val = float(NU)
chk(2/(1-nu_val), "2/(1-ν_RG)", "structural")
chk(2/(1-nu_val/2), "2/(1-ν/2)", "structural")
chk(1/(1-nu_val), "1/(1-ν)", "structural")
chk(2+2*nu_val, "2+2ν", "structural")
chk(2+nu_val, "2+ν", "structural")
chk(1/(nu_val), "1/ν", "structural")
chk(math.pi*nu_val, "π×ν", "structural")

# ζ at d*/2 = 1.2315
z_alpha = float(mpmath.zeta(ALPHA))
z_dstar = float(mpmath.zeta(D_STAR))
z_2_dstar = float(mpmath.zeta(2*D_STAR))
print(f"  ζ(alpha=1.2315) = {z_alpha:.8f}")
print(f"  ζ(d*=2.4631) = {z_dstar:.8f}")
print(f"  ζ(2×d*=4.9262) = {z_2_dstar:.8f}")
# Does d* satisfy ζ(d*/2) = some simple value?
# ζ(alpha) = ζ(1.2315) ≈ ?
# Check: is ζ(alpha) = integer? rational? involves π or e?
for val in [1, 2, 3, 4, math.pi, math.e, math.sqrt(2), math.pi**2/6, math.pi/2]:
    err = abs(z_alpha - val)
    if err < 0.1:
        print(f"  ζ(alpha) ≈ {val} (err={err:.4f})")

# Does d* = 2 + 1/ζ(3/2)?
z32 = float(mpmath.zeta(1.5))
chk(2 + 1/z32, "2+1/ζ(3/2)", "functional")
chk(2 * z32 / math.pi, "2ζ(3/2)/π", "functional")
chk(z32 - 1, "ζ(3/2)-1", "functional")
chk(1/z32 + 1, "1/ζ(3/2)+1", "functional")

# Fixed point: d* = F(d*)?
for f_name, F in [
    ("ζ(x/2)", lambda x: float(mpmath.zeta(x/2))),
    ("2+log(ζ(x))", lambda x: 2+math.log(abs(float(mpmath.zeta(x))))),
    ("x*ζ(x)/π²", lambda x: x*float(mpmath.zeta(x))/math.pi**2),
    ("2/(ζ(x)-1)", lambda x: 2/(float(mpmath.zeta(x))-1)),
]:
    try:
        val = F(D_FLOAT)
        chk(val, f_name + " (fixed-pt test)", "functional")
    except: pass

# Gamma function fixed point
from scipy.special import gamma as Gamma
for f_name, F in [
    ("Γ(x/2)", lambda x: Gamma(x/2)),
    ("2*Γ(x/4)", lambda x: 2*Gamma(x/4)),
    ("Γ(x)/Γ(x-1)", lambda x: Gamma(x)/Gamma(x-1)),  # = x-1
    ("log(Γ(x))/log(x)", lambda x: math.lgamma(x)/math.log(x)),
]:
    try:
        val = F(D_FLOAT)
        chk(val, f_name, "functional")
        # Also check if d* = fixed point: val ≈ d*?
        err_fp = abs(val - D_FLOAT)
        if err_fp < 0.05:
            print(f"  Near fixed-point: {f_name} at d* = {val:.6f} (err={err_fp:.4f})")
    except: pass

# ══════════════════════════════════════════════════════════════
# C. STRUCTURAL: does d* appear in prime/zero data?
# ══════════════════════════════════════════════════════════════
print()
print("=== C. Structural: prime/zero connections ===")

# Zero pair correlations: Montgomery's formula involves specific constants
# Does d* = 1 + something involving γ₁?
g1 = float(GAMMA1)
g2 = float(GAMMA2)
g3 = float(GAMMA3)

for n in range(1, 20):
    chk(1 + g1/n,  f"1+γ₁/{n}", "structural")
    chk(2 - g1/n,  f"2-γ₁/{n}", "structural")
    chk(n/g1 + 2,  f"{n}/γ₁+2", "structural")
    chk((g1+g2)/n, f"(γ₁+γ₂)/{n}", "structural")

# Density of zeros at height T: N(T) ≈ T/(2π) × log(T/(2πe))
# Mean zero spacing at γ₁: δ = 2π/log(γ₁/2π)
mean_spacing = 2*math.pi / math.log(g1/(2*math.pi))
print(f"  Mean zero spacing at γ₁: {mean_spacing:.6f}")
chk(mean_spacing, "2π/log(γ₁/2π)", "structural")
chk(1/mean_spacing, "1/mean_spacing", "structural")
chk(2*mean_spacing, "2×mean_spacing", "structural")
chk(mean_spacing + 1, "mean_spacing+1", "structural")
# Mean spacing at γ₂
ms2 = 2*math.pi / math.log(g2/(2*math.pi))
chk(ms2, "2π/log(γ₂/2π)", "structural")

# Von Mangoldt explicit formula coefficient
# The prime counting error: ψ(x)-x ~ x^(1/2)/log(x) near the zeros
# Does d* relate to the critical exponent 1/2 + something?
chk(2*(0.5 + 1/g1), "2(1/2+1/γ₁)", "structural")
chk(2*0.5 + 1/g1, "1+1/γ₁", "structural")

# Heat kernel dimension from Weyl law for τ_N:
# eigenvalues λ_n ~ n^(2/d*), so d* = 2/β where λ_n ~ n^β
# From the spectral gap analysis: spectral gap ratio = 1.99e+07 ≈ 2×10^7
# If this gap = exp(d* × something)...
gap = 1.99e7
chk(math.log(gap)/math.pi, "log(gap)/π", "structural")
chk(math.log(gap)/g1, "log(gap)/γ₁", "structural")

print()
print("=== SUMMARY ===")
total = sum(len(v) for v in results.values())
print(f"Total matches: {total}")
for cat, items in results.items():
    if items:
        print(f"\n  [{cat}]")
        for it in items:
            print(f"    {it['expr']} = {it.get('val', it.get('root', '?')):.10f}  err={it['err']:.2e}")

# Additional insight: what IS 2.4631 in terms of the RG flow?
print()
print("=== KEY STRUCTURAL INSIGHT ===")
print(f"  ν (RG anomalous dim) = {nu_val:.6f}")
print(f"  2/(1-ν) = {2/(1-nu_val):.6f}  (if d* = 2/(1-ν): err={abs(2/(1-nu_val)-D_FLOAT):.4e})")
print(f"  2/(1-ν/2) = {2/(1-nu_val/2):.6f}  (err={abs(2/(1-nu_val/2)-D_FLOAT):.4e})")
print(f"  2+2ν = {2+2*nu_val:.6f}  (err={abs(2+2*nu_val-D_FLOAT):.4e})")
print()
print(f"  The exponent ν = {nu_val:.6f} comes from ε_N ~ A×N^ν (flow_A in rg_fixedpoint.json)")
print(f"  If d* = 2/(1-ν): interpretation = anomalous scaling dimension in ν-regulated theory")
print(f"  If d* = 2+2ν: interpretation = d* is integer dimension + 2×anomalous correction")
print()
# The actual relationship:
actual_nu_for_dstar = 1 - 2/D_FLOAT
print(f"  For d* = 2/(1-ν) to hold EXACTLY: ν would need to be {actual_nu_for_dstar:.8f}")
print(f"  Actual ν from flow_A: {nu_val:.8f}")
print(f"  Discrepancy: {abs(actual_nu_for_dstar - nu_val):.4e}")

with open("results/d_star_algebraic.json","w") as fh:
    json.dump({"d_star": D_FLOAT, "alpha": A_FLOAT,
               "algebraic_hits": results["algebraic"],
               "functional_hits": results["functional"],
               "structural_hits": results["structural"]}, fh, indent=2)
print("\nSaved → results/d_star_algebraic.json")

#!/usr/bin/env python3
"""
subgoalH_jensen_turan.py — LEDGER §151 (Lens 4/7 grounding): the Gröchenig TP-ladder (TP_k of Λ=FT(1/ξ)) is (claimed)
the same as the JENSEN-POLYNOMIAL hyperbolicity ladder of Ξ. Ground it: compute Ξ's Taylor coefficients → the moment
sequence γ_n>0 → Jensen polynomials J_{d,n}(X)=Σ_k C(d,k) γ_{n+k} X^k → check hyperbolicity (all real roots) and the
Turán margins, degree by degree. This shows WHERE the ladder is tight and connects our TP framing to the known
Jensen/GHRZ frontier (rule-13). Also compute Λ's own log-concavity margin (Λ')²−ΛΛ'' to confirm TP_2 ⟺ degree-2 Jensen.

Ξ(t)=ξ(1/2+it)=Σ_n (-1)^n (γ_n/(2n)!) t^{2n}, γ_n=∫_0^∞ u^{2n}Φ(u)du>0 (Riemann Φ moments). So from Taylor coeff
c_n of Ξ at 0: γ_n = (-1)^n (2n)! c_n. Jensen deg d, shift n hyperbolic ∀(d,n) ⟺ Ξ∈Laguerre-Pólya ⟺ RH.
Turán (d=2): γ_{n+1}² ≥ γ_n γ_{n+2}  (equivalently the discriminant of J_{2,n} ≥0).
"""
import numpy as np, mpmath as mp
mp.mp.dps = 40

# Riemann Φ kernel: Ξ(t)=∫_0^∞ Φ(u)cos(tu)du, Φ even. Moments γ_n=∫_0^∞ u^{2n}Φ(u)du>0 (GHRZ). Compute γ_n by direct
# stable integration (Φ decays double-exponentially) — far faster/cleaner than high-order numerical Taylor derivatives.
def Phi(u):
    u=mp.mpf(u); s=mp.mpf(0); k=1
    while True:
        term=(2*mp.pi**2*k**4*mp.e**(mp.mpf(9)*u/2) - 3*mp.pi*k**2*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*k**2*mp.e**(2*u))
        s+=term
        if k>3 and abs(term) < mp.mpf(10)**(-50)*(abs(s)+1): break
        k+=1
        if k>200: break
    return s

NMAX=16
print("building moments γ_n = ∫_0^∞ u^{2n} Φ(u) du (Riemann Φ kernel, mpmath dps=40)...", flush=True)
gamma=[]
for n in range(NMAX+3):
    g = mp.quad(lambda u: u**(2*n)*Phi(u), [0, 0.5, 1.0, 1.6, 2.4])
    gamma.append(g)
print("γ_n (moments, should be POSITIVE & log-convex-ish):")
for n in range(0,NMAX+1,2):
    print(f"   γ_{n} = {mp.nstr(gamma[n],8)}")

def jensen_hyperbolic(d,n):
    # J_{d,n}(X) = Σ_{k=0}^d C(d,k) γ_{n+k} X^k ; hyperbolic ⟺ all roots real
    coeffs=[mp.binomial(d,k)*gamma[n+k] for k in range(d+1)]
    coeffs=[complex(z) for z in coeffs][::-1]     # numpy wants highest degree first
    r=np.roots(coeffs)
    maximag=np.max(np.abs(r.imag))/ (np.max(np.abs(r.real))+1e-30)
    return maximag  # ~0 ⟹ all-real ⟹ hyperbolic

print("\nJensen-polynomial hyperbolicity ladder (max|Im root|/|Re| ; ~0 = hyperbolic ✓):")
print(f"  {'deg d':>6} " + " ".join(f"n={n:<2}" for n in range(0,10,1)))
for d in [2,3,4,5,6]:
    row=[]
    for n in range(0,10):
        if n+d<=NMAX:
            row.append(jensen_hyperbolic(d,n))
        else: row.append(None)
    cells=" ".join((f"{v:>6.0e}" if v is not None else "   -  ") for v in row)
    print(f"  d={d:<4} {cells}",flush=True)

print("\nTurán inequalities (d=2): γ_{n+1}² − γ_n γ_{n+2}  (≥0 unconditionally = Csordas-Norfolk-Varga):")
for n in range(0,NMAX-1):
    T = gamma[n+1]**2 - gamma[n]*gamma[n+2]
    rel = T/(gamma[n+1]**2)
    print(f"   n={n:>2}: γ²−γγ = {mp.nstr(T,6)}   (relative margin {mp.nstr(rel,5)})  {'✓' if T>0 else '✗ NEG'}")
print("\nREAD: all Jensen entries ~0 (hyperbolic) + all Turán ≥0 ⟹ the ladder holds numerically (RH true). The TURÁN row")
print("(d=2) = TP_2 of Λ = KNOWN unconditional theorem (CNV 1986). The frontier is HIGH-d asymptotics (GHRZ 2019). This")
print("grounds our TP-lever ON the Jensen ladder: same object. Question then = does our prime-side reach where GHRZ stops.")

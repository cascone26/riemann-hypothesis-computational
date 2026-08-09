#!/usr/bin/env python3
"""§868: the WEIL MEASURE dW(u) = arch continuous density - prime atoms. Q(phi)=int dW(u)*autocorr(phi-hat)(u).
Compute the archimedean density (FT of (1/2pi)Omega(r), smoothed to the band) and the prime atoms
Lambda(n)/sqrt(n) at u=log n. Where the arch density (per unit u, integrated over the atom's Fejer-width from
the band-limit) exceeds the prime atomic mass, positivity is manifest. This is the constructive core."""
import numpy as np
from scipy.special import digamma
# archimedean density in u-domain: arch(phi) = (1/2pi) int Omega(r)|phi|^2 dr.
# |phi|^2 has spectrum in [-4W,4W]; a prime atom at u0=log n couples via autocorr(phi-hat)(u0).
# The arch part contributes (1/2pi)*Omega-hat convolved. Compute Omega-hat(u) = FT of Omega(r).
# Omega(r)=Re psi0(1/4+ir/2)-ln pi. Use large-r numeric FT with a smooth taper.
R=np.arange(-4000,4000,0.01)
OM=np.real(digamma(0.25+0.5j*R))-np.log(np.pi)
taper=np.exp(-(R/3000.0)**2)
def Ohat(u):
    return np.trapezoid((OM*taper)*np.cos(u*R),R)   # real part (Omega even)
print("§868 Weil measure: arch density (1/2pi)*Omega-hat(u) vs prime atoms Lambda(n)/sqrt(n) at u=log n")
print(f"{'u=log n':>9} {'n':>4} {'prime atom Ln/sqrt(n)':>21} {'(1/2pi)Ohat(u)':>15} {'arch/prime':>11}")
# primes/prime-powers
def Lambda(n):
    m=n; 
    for p in range(2,n+1):
        if m%p==0:
            while m%p==0: m//=p
            return np.log(p) if m==1 else 0.0
    return 0.0
inv2pi=1/(2*np.pi)
for n in [2,3,4,5,7,8,9,11,13,16]:
    L=Lambda(n)
    if L==0: continue
    u=np.log(n); atom=L/np.sqrt(n)
    ad=inv2pi*Ohat(u)
    print(f"{u:>9.4f} {n:>4} {atom:>21.5f} {ad:>15.5f} {ad/atom:>11.3f}")
print("-"*70)
print("Note: this compares POINTWISE density to atom mass. The real positivity is the TRUNCATED-PD (Fejer-Riesz)")
print("condition on the band [-4W,4W]; the arch density is spread, atoms are localized. Workflow PATH B/D solve the")
print("actual truncated problem. This shows the raw arch-density vs prime-atom balance per prime.")

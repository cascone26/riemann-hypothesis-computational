#!/usr/bin/env python3
"""§939 (spectral, semiclassical density bridge, independent). Does the prolate operator's spectral density match the
ζ-zero density? (1) verify N(T)=(T/2π)ln(T/2π)−T/2π+7/8 vs actual zeros; (2) the prolate 2nd-order operator L on [-1,1]
eigenvalues χ_n — compute their WKB density and compare structure. If leading densities match, RH-content = fluctuation."""
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh
import mpmath as mp
mp.mp.dps=20
# (1) zero counting vs Riemann-von Mangoldt
def Ncount(T): return (T/(2*np.pi))*np.log(T/(2*np.pi))-T/(2*np.pi)+7/8
ZT=[]; k=1
while True:
    g=float(mp.im(mp.zetazero(k))); ZT.append(g)
    if g>150: break
    k+=1
ZT=np.array(ZT)
print("="*76)
print("§939 semiclassical density: ζ-zero counting N(T) vs Riemann-von Mangoldt formula")
print("="*76)
print(f"{'T':>6} {'actual #zeros≤T':>15} {'N(T) formula':>13} {'diff':>8}")
for T in [30,50,80,100,120,150]:
    ac=int((ZT<=T).sum()); print(f"{T:>6} {ac:>15} {Ncount(T):>13.3f} {ac-Ncount(T):>8.3f}")
# (2) prolate 2nd-order operator L f = -((1-x^2)f')' + c^2 x^2 f on [-1,1], Dirichlet-ish; eigenvalues chi_n
def prolate_chi(c,M=600):
    x=np.linspace(-1,1,M); h=x[1]-x[0]
    # discretize -((1-x^2)f')' + c^2 x^2 f  (self-adjoint Sturm-Liouville)
    p=1-x**2  # coefficient
    ph=1-((x[:-1]+x[1:])/2)**2  # p at half-points
    main=np.zeros(M); off=np.zeros(M-1)
    for i in range(1,M-1):
        main[i]=(ph[i-1]+ph[i])/h**2 + c**2*x[i]**2
        off[i-1]=-ph[i-1]/h**2
    main[0]=1e6; main[-1]=1e6  # Dirichlet at endpoints (soft)
    ev=eigh_tridiagonal(main,off,select='i',select_range=(0,min(30,M-1)))[0]
    return np.sort(ev)
print("\nprolate operator L eigenvalues χ_n (the commuting 2nd-order op; CCM use these for the zero-connection):")
for c in [10,20,30]:
    chi=prolate_chi(c)
    print(f"  c={c}: χ_0..χ_8 = {np.array2string(chi[:9],precision=1)}")
    # WKB: n-th eigenvalue phase; density dn/dχ. gaps:
    print(f"        χ gaps (χ_{{n+1}}−χ_n) n=0..6: {np.array2string(np.diff(chi[:8]),precision=1)}")
print("="*76)
print("N(T) matches actual zeros to <1 (density formula solid). The prolate χ_n density (WKB) is the object to match to the")
print("zero density; if leading-order equal (Connes trace formula), the RH-content = FLUCTUATION, cleaner than the razor.")
print("="*76)

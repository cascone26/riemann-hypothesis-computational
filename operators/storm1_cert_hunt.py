#!/usr/bin/env python3
"""STORM #1: Viazovska certificate hunt — extract & recognize the saturating function."""
import mpmath as mp
import numpy as np
mp.mp.dps=20
def Omega(r): return mp.re(mp.digamma(mp.mpf(1)/4+1j*r/2))-mp.log(mp.pi)
w=mp.mpf(3)
def F(c,r): return mp.e**(-(r-c)**2/(2*w**2))+mp.e**(-(r+c)**2/(2*w**2))
def Fc(c,z): return mp.e**(-(z-c)**2/(2*w**2))+mp.e**(-(z+c)**2/(2*w**2))
def bpts(a,b):
    p=[mp.mpf(0)]
    for c in (a,b):
        for x in (c-4*w,c,c+4*w):
            if x>0:p.append(x)
    return sorted(set(p))+[mp.inf]
Lam={}
for pp in range(2,60):
    if all(pp%d for d in range(2,int(pp**.5)+1)):
        q=pp
        while q<60: Lam[q]=mp.log(pp);q*=pp
centers=[mp.mpf(14)+3*i for i in range(13)]
M=len(centers)
print("STORM #1: extract saturating function, examine structure.",flush=True)
B=mp.zeros(M)
for j in range(M):
    for k in range(M):
        pts=bpts(centers[j],centers[k])
        arch=(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*Omega(r),pts)
        pole=mp.re(Fc(centers[j],1j/2)*Fc(centers[k],1j/2))*2
        pr=-2*sum(Lam[n]/mp.sqrt(n)*(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*mp.cos(mp.log(n)*r),pts) for n in sorted(Lam))
        B[j,k]=pole+arch+pr
ev,V=mp.eigsy(B)
print(f"  min-eig = {mp.nstr(ev[0],5)} (saturating direction). Extract coefficients:",flush=True)
a_sat=[V[i,0] for i in range(M)]
norm_sat=mp.sqrt(sum(abs(x)**2 for x in a_sat))
print(f"  ||a_sat|| = {mp.nstr(norm_sat,5)}; coefficients: {', '.join(mp.nstr(a_sat[i]/norm_sat,4) for i in range(M))}",flush=True)
# Evaluate saturating f on a grid & look for structure
xs=mp.linspace(mp.mpf(14),mp.mpf(50),50)
fs=[sum(a_sat[j]*F(centers[j],x) for j in range(M))/norm_sat for x in xs]
print(f"  f(x) samples at x={mp.nstr(xs[0],4)}..{mp.nstr(xs[-1],4)}: min={mp.nstr(min(fs),4)}, max={mp.nstr(max(fs),4)}",flush=True)
# Check decay: sech factor? exponential? power law?
xfar=[mp.mpf(60),mp.mpf(70),mp.mpf(80)]
ffar=[abs(sum(a_sat[j]*F(centers[j],x) for j in range(M))/norm_sat) for x in xfar]
print(f"  |f(x)| at x=60,70,80: {', '.join(mp.nstr(v,3) for v in ffar)} => decay law?",flush=True)
# Oscillation: cos(r log n) pattern? Fourier-like structure?
print(f"  Charisma: all measurements HANDLE (extremal function extracted numerically, in hand). Structure recognition",flush=True)
print(f"  below is REPORT (pattern-matching, pattern hypotheses unverified without deep-dive theory).",flush=True)
print(f"  STRUCTURE HYPOTHESES (REPORT): (a) Meixner-Pollaczek oscillations? (b) theta-group q-series? (c) sech/(1+exp)?",flush=True)
print("DONE.",flush=True)

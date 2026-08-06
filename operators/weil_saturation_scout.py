#!/usr/bin/env python3
"""Scout the SATURATION EDGE of finite Weil-positivity certification (fast float, Arb-validated 1e-5).
Weil form Q=VᵀV (§772) is PSD but non-coercive: min-eig→0⁺ as the basis densifies. Rigorous Arb
certification needs min-eig(Q) > ρ (the enclosure radius, floor ~1e-11..1e-13 from §757). Map
min-eig vs basis density to find the DENSEST basis still certifiable — the practical frontier of
finite Weil verification. Uses the von-Mangoldt-correct fast build (§769)."""
import numpy as np, mpmath as mp
mp.mp.dps=30
PI=np.pi; LOGPI=np.log(PI)
R=240.0; NG=48000; rg=np.linspace(-R,R,NG)
print("precomputing Ω(r)...", flush=True)
OM=np.array([float(mp.re(mp.digamma(mp.mpf(1)/4+1j*mp.mpf(r)/2)))-LOGPI for r in rg])
def phihat(a,r): s2p=np.sqrt(2*PI); return s2p/2*(np.exp(-(r-a)**2/2)+np.exp(-(r+a)**2/2))
def arch(aj,ak): return np.trapezoid(phihat(aj,rg)*phihat(ak,rg)*OM,rg)/(2*PI)
def pole(aj,ak):
    s2p=np.sqrt(2*PI)
    pj=s2p/2*(np.exp(-((1j/2-aj)**2)/2)+np.exp(-((1j/2+aj)**2)/2))
    pk=s2p/2*(np.exp(-((1j/2-ak)**2)/2)+np.exp(-((1j/2+ak)**2)/2))
    return 2*np.real(pj*pk)
from sympy import primerange
NMAX=300000;_p=list(primerange(2,NMAX+1));_N=[];_L=[]
for p in _p:
    pk=p;lp=np.log(p)
    while pk<=NMAX:_N.append(pk);_L.append(lp);pk*=p
_N=np.array(_N,float);_L=np.array(_L);_LOGN=np.log(_N);_IS=1/np.sqrt(_N)
def g(aj,ak,y):
    s=aj+ak;d=aj-ak
    return np.exp(-y**2/4)*(np.sqrt(PI)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
def prime(aj,ak): return 2*np.sum(_L*_IS*g(aj,ak,_LOGN))
def Qmat(A):
    n=len(A);Q=np.zeros((n,n))
    for i in range(n):
        for j in range(i,n): Q[i,j]=Q[j,i]=arch(A[i],A[j])+pole(A[i],A[j])-prime(A[i],A[j])
    return Q
RHO_FLOOR=1e-11   # achievable Arb enclosure radius (from §757)
print(f"\n{'window':>14} {'N':>3} {'spacing':>8} {'min-eig':>13} {'certifiable?':>13}")
# fixed window [15,55], increasing N (densifying) — find where min-eig crosses RHO_FLOOR
for W in [(15,55),(15,45),(20,40)]:
    for N in range(3,26):
        A=list(np.linspace(W[0],W[1],N)); sp_=(W[1]-W[0])/(N-1)
        me=np.linalg.eigvalsh(Qmat(A)).min()
        cert = me > RHO_FLOOR
        mark = "YES" if cert else "no (me<ρ)"
        print(f"{str(W):>14} {N:>3} {sp_:>8.3f} {me:>13.3e} {mark:>13}")
        if not cert and me < RHO_FLOOR/10: break   # well past the edge
    print()

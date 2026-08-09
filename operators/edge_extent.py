#!/usr/bin/env python3
"""§911: CRITICAL — does the corner (edge) region stay BOUNDED as 4W→∞? For each 4W, sweep windows climbing in height
and find the MAX height at which a corner (#neg(M_σ)=1) still appears near capacity. If max-corner-height stays bounded
or grows slowly (log) → the edge region is ~finite → finite-proof path (§910). If it grows fast → corner not confined."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3200; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def hasneg(a,b,fourW,frac,dr=0.006):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+50),(b+50),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    gw=np.linalg.eigvalsh(Gm)
    if gw[0]<=1e-13: return None
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    gw2,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw2))@gv.T
    ev=eigh(B@(arch-prime)@B,eigvals_only=True)
    return int((ev<-1e-10).sum()), float(ev[0])
print("="*90)
print("§911 does the corner region stay BOUNDED as 4W grows? MAX corner-height per 4W (frac 0.85)")
print("="*90)
print(f"{'4W':>5} {'max corner height':>18} {'(deepest window, μ)':>28}")
for fourW in [4.0,5.0,6.0,7.0,8.0]:
    maxh=0; deepest=None
    for (a,b) in [(14,30),(14,44),(30,50),(30,60),(45,70),(60,90),(80,115),(100,140),(130,175),(170,220)]:
        r=hasneg(a,b,fourW,0.85)
        if r is None: continue
        nn,ev0=r
        if nn>=1:
            h=(a+b)/2
            if h>maxh: maxh=h; deepest=(a,b,ev0)
    ds=f"[{deepest[0]},{deepest[1]}] μ={-deepest[2]:.1e}" if deepest else "none"
    print(f"{fourW:>5.1f} {maxh:>18.0f} {ds:>28}")
print("="*90)
print("If max corner-height ~const or ~log(4W) → edge region ~finite → finite-proof path. If ~linear in 4W → grows.")
print("="*90)

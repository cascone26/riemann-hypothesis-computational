#!/usr/bin/env python3
"""§904 (independent): does the corner #neg(M_σ)=1 EXIST for small 4W? Sweep 4W from 1→7, and for each, scan
windows/heights near capacity, reporting the MAX #neg(M_σ) found and the min margin. If #neg(M_σ)=0 for ALL configs
below some 4W threshold → a provable sub-threshold regime (M_σ⪰0, positivity trivial). Cross-checks workflow Angle 1."""
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
NMAX=6000; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def negcount(a,b,fourW,frac,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
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
    return int((ev<-1e-11).sum()), float(ev[0])
print("="*92)
print("§904 does the corner #neg(M_σ)=1 exist for small 4W? (2π≈6.283). MAX #neg over windows×heights×frac per 4W")
print("="*92)
print(f"{'4W':>5} {'maxneg':>7} {'min(minM_σ)':>13}  configs with #neg>0 (window,height)")
for fourW in [1.0,2.0,3.0,4.0,5.0,6.0,6.3,7.0]:
    maxneg=0; minev=np.inf; hits=[]
    for (a,b) in [(14,44),(14,34),(20,50),(30,60),(45,75),(8,24),(10,26)]:
        for frac in [0.7,0.85,0.95]:
            r=negcount(a,b,fourW,frac)
            if r is None: continue
            nn,ev0=r; maxneg=max(maxneg,nn); minev=min(minev,ev0)
            if nn>0: hits.append(f"[{a},{b}]")
    hs=",".join(sorted(set(hits))[:6])
    print(f"{fourW:>5.1f} {maxneg:>7} {minev:>13.3e}  {hs}")
print("="*92)
print("If maxneg=0 for all 4W below a threshold → provable sub-threshold (M_σ⪰0, no corner). Watch vs 2π≈6.28.")
print("="*92)

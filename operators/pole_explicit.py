#!/usr/bin/env python3
"""§907: characterize the POLE COVERAGE π_w=p(v·w)² explicitly (the tractable half of RH-corner ⟺ π_w≥μ). The pole
= 2Re(f fᴴ), f_j = Fejér_j(0.5i±c) = test-function value at the pole point s=1 (r=i/2) — EXPLICIT archimedean data.
Compute p (pole strength), the pole direction v, and how (v·w)² and π_w scale with window height/bandwidth. Goal: is
π_w bounded below by an explicit closed form? (μ is the RH-hard half; π_w is not.)"""
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
def run(a,b,fourW,frac=0.7,dr=0.003):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+70),(b+70),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    Msw=B@(arch-prime)@B; polew=B@pole@B
    ev,V=eigh(Msw); mu=-ev[0]; w=V[:,0]
    pev,pvec=eigh(polew); p=pev[-1]; v=pvec[:,-1]
    # ||fp|| (raw pole vector norm before whitening) and its decay
    fpnorm=np.linalg.norm(fp)
    return N,Ncap,mu,p,(v@w)**2,p*(v@w)**2,fpnorm,p*(v@w)**2/mu
print("="*104)
print("§907 pole coverage π_w=p(v·w)² — the EXPLICIT (archimedean) half of RH-corner ⟺ π_w≥μ. scaling vs height/4W")
print("="*104)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'μ':>10} {'p(pole)':>10} {'(v·w)²':>8} {'π_w':>10} {'‖f(i/2)‖':>9} {'π_w/μ':>7}")
for (a,b,fourW) in [(14,44,4.0),(14,34,4.0),(20,50,4.0),(30,60,4.0),(14,44,3.0),(14,44,5.0),(45,75,4.0)]:
    N,Ncap,mu,p,vw2,pw,fpn,ratio=run(a,b,fourW)
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {mu:>10.3e} {p:>10.3e} {vw2:>8.4f} {pw:>10.3e} {fpn:>9.2e} {ratio:>7.3f}")
print("="*104)
print("p (pole strength) is EXPLICIT archimedean data (Fejér at r=i/2). If π_w has a clean lower bound vs μ → half a proof.")
print("Note ‖f(i/2)‖ = raw pole-vector norm; p = its whitened strength. Both zero-free/explicit.")
print("="*104)

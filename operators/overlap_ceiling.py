#!/usr/bin/env python3
"""§903: the finest RH-safety test. Razor ≈ g = 1 − (w·a_dom)² (§901), w=M_σ neg eigvec, a_dom=dominant near-null
A-mode. Does (w·a_dom)² stay BOUNDED BELOW 1 (RH-safe, a real ceiling) or → 1 (razor→0) as the corner deepens
(frac→1, higher 4W)? Compute (w·a_dom)² and the exact g across corners at frac→0.97, track cond(G). Also report the
'ceiling gap' 1−(w·a_dom)². This is §900's floor test at the FINEST (true-margin) resolution."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh, pinvh
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=6000; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def build(a,b,fourW,N,dr=0.003):
    W=fourW/4; RG=np.arange(-(b+70),(b+70),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr; cond=np.linalg.cond(Gm)
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    return B@(arch-prime)@B, B@pole@B, cond
print("="*104)
print("§903 does (w·a_dom)² stay < 1 (RH ceiling) or → 1 (razor→0)? finest-resolution safety test, frac→0.97")
print("="*104)
for (a,b,fourW) in [(14,44,4.0),(14,44,5.0),(14,34,4.0),(20,50,4.0)]:
    Ncap=(b-a)*2*(fourW/4)/np.pi
    print(f"\n[{a},{b}] 4W={fourW} Ncap={Ncap:.1f}:  {'N':>3} {'frac':>5} {'cond(G)':>8} {'(w·a_dom)²':>11} {'1−(w·a_dom)²':>12} {'exact g':>10} {'margin':>11}")
    for frac in [0.55,0.7,0.85,0.94,0.97]:
        N=max(4,int(round(frac*Ncap)))
        Msw,polew,cond=build(a,b,fourW,N)
        ev,V=eigh(Msw); nn=int((ev<0).sum())
        if nn<1: continue
        mu=-ev[0]; w=V[:,0]; Mplus=V[:,1:]@np.diag(ev[1:])@V[:,1:].T
        A=Mplus+polew; al,ak=eigh(A)
        contrib=np.array([(w@ak[:,k])**2/al[k] for k in range(N)]); kdom=np.argmax(contrib)
        wa2=(w@ak[:,kdom])**2
        g=1-mu*contrib.sum(); margin=float(eigh(Msw+polew,eigvals_only=True)[0])
        print(f"{'':>{len(f'[{a},{b}] 4W={fourW} Ncap={Ncap:.1f}:')}}  {N:>3} {frac:>5.2f} {cond:>8.1e} {wa2:>11.6f} {1-wa2:>12.3e} {g:>+10.3e} {margin:>+11.3e}")
print("="*104)
print("(w·a_dom)² bounded below 1 with a stable gap as frac→1 → RH-safe ceiling (proof target). →1 → razor asymptotically zero.")
print("="*104)

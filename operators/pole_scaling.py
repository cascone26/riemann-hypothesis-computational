#!/usr/bin/env python3
"""§915: POLE-SCALING sensitivity α* (new lens, dual to §856 prime-scaling). Scale the pole: Q(α)=arch+α·pole−prime.
α=1 is the real Weil form. α* = smallest α making Q(α)⪰0. α*<1 → the real pole OVER-covers (slack = 1−α*); α*=1 →
razor-exact; α*>1 → would need MORE pole (RH would fail). Map α* across the bottom razor: how much pole slack is there?"""
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
def astar(a,b,fourW,frac=0.85,dr=0.005):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+50),(b+50),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    if np.linalg.eigvalsh(Gm)[0]<=1e-13: return None
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    AP=arch-prime
    def minev(al): return float(eigh(AP+al*pole,Gm,eigvals_only=True)[0])
    if minev(0)>=0: return (N,0.0)  # no pole needed (M_σ⪰0)
    lo,hi=0.0,1.0
    while minev(hi)<0 and hi<50: lo,hi=hi,hi*2
    for _ in range(40):
        m=(lo+hi)/2
        if minev(m)>=0: hi=m
        else: lo=m
    return (N,hi)
print("="*88)
print("§915 pole-scaling α* (smallest α with arch+α·pole−prime⪰0). α*=1 razor-exact; α*<1 pole over-covers (slack=1−α*)")
print("="*88)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'α*':>8} {'pole slack 1−α*':>16}")
for fourW in [4.0,5.0,6.0,7.0]:
    for (a,b) in [(14,30),(14,44),(20,44),(24,48)]:
        r=astar(a,b,fourW)
        if r is None: continue
        N,al=r
        if al==0.0: print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {'--':>8} {'M_σ⪰0 (no pole)':>16}"); continue
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {al:>8.4f} {1-al:>16.4f}")
print("="*88)
print("α* just below 1 (slack small +) = real pole barely covers (razor). α* ≪ 1 = comfortable. α* > 1 = would fail (never, if RH).")
print("="*88)

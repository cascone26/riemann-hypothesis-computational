#!/usr/bin/env python3
"""§909: how LOCAL is the adjacent-zero spillover (§908)? For the corner worst-direction w, compute the razor
contribution 2w(γ)² per zero, ordered by γ, and the cumulative fraction. How many adjacent above-zeros account for
90% of the margin? Does contribution DECAY with distance above the window edge b? If localized to nearest ~5-10 zeros
with a clean decay → a local extremal object (Beurling–Selberg boundable)."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import mpmath as mp
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=6000; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
mp.mp.dps=25
ZT=[]; k=1
while True:
    g=float(mp.im(mp.zetazero(k))); ZT.append(g)
    if g>200: break
    k+=1
ZT=np.array(ZT)
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
    ev,V=eigh(arch+pole-prime,Gm); coef=V[:,0]
    def wval(r): return sum(coef[j]*(fej(r-CEN[j])+fej(r+CEN[j])) for j in range(N))
    contrib=np.array([2*wval(g)**2 for g in ZT])
    return contrib
print("="*100)
print("§909 how local is the spillover? razor contribution 2w(γ)² per zero, by distance from window edge b")
print("="*100)
for (a,b,fourW) in [(14,44,4.0),(20,50,4.0),(14,34,4.0)]:
    c=run(a,b,fourW); tot=c.sum()
    below=c[ZT<a].sum(); inwin=c[(ZT>=a)&(ZT<=b)].sum(); above=c[ZT>b].sum()
    print(f"\n[{a},{b}] 4W={fourW}: total 2Σw(γ)²={tot:.3e}")
    print(f"  below window ({(ZT<a).sum()} zeros): {100*below/tot:.1f}%   in-window ({((ZT>=a)&(ZT<=b)).sum()}): {100*inwin/tot:.1f}%   above ({(ZT>b).sum()}): {100*above/tot:.1f}%")
    # cumulative over above-zeros nearest-first
    az=ZT[ZT>b]; ac=c[ZT>b]; cum=np.cumsum(ac)/tot
    print(f"  nearest above-zeros γ (dist from b): "+", ".join(f"{g:.1f}(+{g-b:.0f}):{100*ac[i]/tot:.0f}%" for i,g in enumerate(az[:8])))
    n90=np.searchsorted(cum,0.9*(above/tot))+1
    print(f"  cumulative above-% reaches 90%-of-above by the nearest {n90} above-zeros")
print("="*100)
print("If contribution decays with distance above b and ~5-10 nearest above-zeros give ~90% → LOCAL extremal object.")
print("below% vs above%: corner windows are at the bottom edge of the zero spectrum (few/no zeros below).")
print("="*100)

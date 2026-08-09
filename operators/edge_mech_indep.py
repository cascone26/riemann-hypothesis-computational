#!/usr/bin/env python3
"""§910 (independent): is the corner an EDGE-of-spectrum effect? Sweep windows climbing from the bottom (h~19) into the
interior (h~200) at fixed 4W=4, near-capacity. Track #neg(M_σ), corner depth μ, margin, and the spillover split
(below/in/above window). Hypothesis: bottom-edge windows have DEEP one-sided (upward) spillover corners; interior
windows have balanced two-sided spillover → shallower/no corner. Cross-checks workflow Angle D."""
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
    if g>400: break
    k+=1
ZT=np.array(ZT)
def run(a,b,fourW,frac,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+80),(b+80),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
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
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw2,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw2))@gv.T
    ev=eigh(B@(arch-prime)@B,eigvals_only=True); nneg=int((ev<-1e-11).sum()); mu=-ev[0] if nneg>=1 else None
    # spillover of worst direction
    evQ,VQ=eigh(arch+pole-prime,Gm); coef=VQ[:,0]; margin=evQ[0]
    def wval(r): return sum(coef[j]*(fej(r-CEN[j])+fej(r+CEN[j])) for j in range(N))
    c2=np.array([2*wval(g)**2 for g in ZT]); tot=c2.sum()
    below=c2[ZT<a].sum()/tot; above=c2[ZT>b].sum()/tot
    return N,nneg,mu,margin,below*100,above*100,(ZT<a).sum()
print("="*100)
print("§910 corner as edge effect? windows climbing bottom→interior (4W=4, frac 0.85). μ, margin, spillover below/above")
print("="*100)
print(f"{'window':>12} {'ht':>4} {'N':>3} {'#neg':>4} {'μ (depth)':>11} {'margin':>11} {'below%':>7} {'above%':>7} {'#z<a':>5}")
for (a,b) in [(14,44),(30,60),(60,90),(100,130),(150,180),(200,230)]:
    r=run(a,b,4.0,0.85)
    if r is None: continue
    N,nn,mu,margin,bel,abv,zb=r; h=(a+b)//2
    mus=f"{mu:.3e}" if mu else "  (PSD)  "
    print(f"{f'[{a},{b}]':>12} {h:>4} {N:>3} {nn:>4} {mus:>11} {margin:>+11.3e} {bel:>7.1f} {abv:>7.1f} {zb:>5}")
print("="*100)
print("If bottom-edge: deep μ + ~100% above-spillover; interior: shallower/PSD + balanced below≈above → EDGE mechanism.")
print("="*100)

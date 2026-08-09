#!/usr/bin/env python3
"""§908: the DUAL/zeros view. Explicit formula: Q(φ)=2Σ_γ φ(γ)² (sum over nontrivial-zero heights γ). The razor
margin = 2Σ_γ w(γ)² for the worst direction w — a band-limited EXTREMAL problem over the zeros (Beurling–Selberg /
Carneiro–Littmann–Vaaler). Verify the identity margin=2Σw(γ)²/‖w‖²_G against ACTUAL zeros, and characterize which
zeros dominate the razor & how w behaves on the zero set (does it nearly-annihilate them?)."""
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
# zeros up to height ~120
ZT=[]; k=1
while True:
    g=float(mp.im(mp.zetazero(k))); ZT.append(g)
    if g>120: break
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
    Q=arch+pole-prime
    ev,V=eigh(Q,Gm); margin=ev[0]; coef=V[:,0]   # worst direction w (gen-eigvec), Q(w)=margin, ‖w‖_G=1
    # evaluate w(r)=Σ coef_j (fej(r-c)+fej(r+c)) at the zeros (use both ±γ symmetric; the Fejér basis already even)
    def wval(r):
        return sum(coef[j]*(fej(r-CEN[j])+fej(r+CEN[j])) for j in range(N))
    zvals=np.array([wval(g) for g in ZT])
    zerosum=2*np.sum(zvals**2)   # should equal Q(w)=margin (‖w‖_G=1)
    # dominant zeros
    contrib=2*zvals**2; order=np.argsort(contrib)[::-1]
    return N,margin,zerosum,ZT[order[:5]],contrib[order[:5]],(ZT>=a)&(ZT<=b)
print("="*100)
print("§908 dual/zeros view: margin = 2Σ_γ w(γ)²? which zeros dominate the razor? (Beurling–Selberg extremal frame)")
print("="*100)
for (a,b,fourW) in [(14,44,4.0),(14,34,4.0),(20,50,4.0)]:
    N,margin,zsum,topz,topc,inwin=run(a,b,fourW)
    nz=int(((ZT>=a)&(ZT<=b)).sum())
    print(f"\n[{a},{b}] 4W={fourW} N={N}: #zeros in window={nz}")
    print(f"  margin(Q,gen-eig) = {margin:+.4e}")
    print(f"  2Σ_γ w(γ)² (actual zeros to h120) = {zsum:+.4e}   ratio to margin = {zsum/margin:.3f}")
    print(f"  top-5 razor-dominating zeros γ: {np.array2string(topz,precision=2)}  contributions: {np.array2string(topc,precision=2)}")
print("="*100)
print("If 2Σw(γ)² ≈ margin → identity confirmed; the razor is which zeros w can't annihilate. Dominant γ near window edges?")
print("="*100)

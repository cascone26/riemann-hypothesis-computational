#!/usr/bin/env python3
"""§862: Szego connection — does the FULL Weil eigenvalue distribution match the symbol S(r)=(1/2pi)Omega-P
distribution over the window? Szego's thm: eigenvalues of a compressed multiplication operator distribute
like the symbol as N->inf. Compare the Q generalized-eigenvalue histogram to the symbol histogram over the
window, at below-capacity and near-capacity N."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def analyze(a,b,W,N,dr=0.004):
    fourW=4*W; RGRID=np.arange(-(b+120),(b+120),dr); OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
    CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    w=eigh(pole+arch-prime,G,eigvals_only=True)
    # symbol S(r) over window
    Psym=np.zeros_like(RGRID)
    for u,ww in zip(ln,lw): Psym+=ww*np.cos(u*RGRID)
    Psym*=(1/np.pi); S=OMEGA/(2*np.pi)-Psym; Sw=S[(RGRID>=a)&(RGRID<=b)]
    return w, Sw
print("="*88); print("§862 Szego: Weil eigenvalue distribution vs symbol S(r) distribution over the window"); print("="*88)
for (a,b,fourW,N) in [(50,80,3.0,8),(50,80,3.0,13),(50,80,3.0,18)]:
    W=fourW/4; w,Sw=analyze(a,b,W,N); Ncap=(b-a)*2*W/np.pi
    print(f"\n[{a},{b}] 4W={fourW} N={N} (N/cap={N/Ncap:.2f}):")
    print(f"  eig:    min {w.min():+.4f}  med {np.median(w):+.4f}  max {w.max():+.4f}  mean {w.mean():+.4f}  #neg {int((w<0).sum())}")
    print(f"  symbol: min {Sw.min():+.4f}  med {np.median(Sw):+.4f}  max {Sw.max():+.4f}  mean {Sw.mean():+.4f}  frac<0 {(Sw<0).mean():.2f}")
    print(f"  eig quartiles:    {np.percentile(w,[10,25,50,75,90]).round(4)}")
    print(f"  symbol quartiles: {np.percentile(Sw,[10,25,50,75,90]).round(4)}")
print("="*88)
print("If eig distribution -> symbol distribution as N->cap (esp. eig min approaching symbol min, negatives appearing)")
print("=> Szego holds; the band-limiting (below cap) is exactly what keeps eigs positive despite negative symbol.")

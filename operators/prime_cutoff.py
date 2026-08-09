#!/usr/bin/env python3
"""§917: how does the razor BUILD as primes are added? At [14,44] 4W=4, compute λmin(arch−prime) and full margin as the
prime sum is truncated at n≤K (K=2,3,5,10,20,54=full). Which primes create the negative eigenvalue? Tests whether the
bottom razor is controlled by a finite low-prime set (finite-argument relevance)."""
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
def build(a,b,fourW,frac=0.85,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    return P,Gm,arch,pole,RG,OM,N,fourW
def primemat(P,RG,K,fourW,dr=0.004):
    N=P.shape[0]; mask=(LOGN<fourW)&(NS<=K); ln=LOGN[mask]; lw=LAMW[mask]
    if len(ln)==0: return np.zeros((N,N))
    cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    return 0.5*(prime+prime.T)
print("="*80)
print("§917 how the razor builds: λmin(arch−prime) & full margin vs prime cutoff K, [14,44] 4W=4")
print("="*80)
P,Gm,arch,pole,RG,OM,N,fourW=build(14,44,4.0)
print(f"N={N}, full prime cutoff n<e^4≈{int(np.exp(4))}")
print(f"{'K (primes n≤K)':>16} {'#terms':>7} {'λmin(arch−prime)':>17} {'full margin':>13}")
for K in [1,2,3,4,5,7,10,15,20,30,54,200]:
    prime=primemat(P,RG,K,fourW)
    lmin=float(eigh(arch-prime,Gm,eigvals_only=True)[0])
    marg=float(eigh(arch+pole-prime,Gm,eigvals_only=True)[0])
    nt=int(((LOGN<fourW)&(NS<=K)).sum())
    print(f"{K:>16} {nt:>7} {lmin:>+17.4e} {marg:>+13.4e}")
print("="*80)
print("λmin crosses negative at the K where enough primes accumulate → those primes create the razor.")
print("="*80)

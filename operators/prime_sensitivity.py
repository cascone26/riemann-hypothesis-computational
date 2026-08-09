#!/usr/bin/env python3
"""§916: which PRIMES control the bottom razor? On the worst direction w, prime(w)=Σ_{log n<4W}(2Λ(n)/√n)·G_n(w),
G_n(w)=(1/2π)∫w(r)²cos(r log n)dr. Decompose prime(w) per prime-power n; which few n dominate the defect μ=prime(w)−arch(w)?
Reveals the arithmetic content of the razor (which primes push arch−prime negative at the bottom edge)."""
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
def run(a,b,fourW,frac=0.85,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; nn=NS[mask]
    cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    ev,V=eigh(arch-prime,Gm); coef=V[:,0]/np.sqrt(V[:,0]@Gm@V[:,0]); mu=-ev[0]
    wr=(coef[:,None]*P).sum(axis=0)  # w(r)
    archw=(wr**2*OM).sum()*dr/(2*np.pi)
    # per-prime contribution to prime(w) = 2Λ(n)/√n · (1/2π)∫ w² cos(r log n) dr
    Gn=np.array([(wr**2*np.cos(RG*l)).sum()*dr/(2*np.pi) for l in ln])
    contrib=2*lw*Gn  # per-n contribution to prime(w)
    return N,mu,archw,contrib.sum(),nn,ln,contrib
print("="*90)
print("§916 which primes control the bottom razor? per-prime-power contribution to prime(w); μ=prime(w)−arch(w)")
print("="*90)
for (a,b,fourW) in [(14,44,4.0),(14,30,4.0),(14,44,5.0)]:
    N,mu,archw,primew,nn,ln,contrib=run(a,b,fourW)
    print(f"\n[{a},{b}] 4W={fourW} N={N}: arch(w)={archw:.4f} prime(w)={primew:.4f} μ=prime−arch={primew-archw:.3e}")
    order=np.argsort(np.abs(contrib))[::-1]
    print("  top prime-powers by |contribution to prime(w)|:")
    for i in order[:8]:
        print(f"    n={nn[i]:3d} (log n={ln[i]:.3f}) Λ/√n={0.5*2*np.exp(-0)*0:.0f}  contrib={contrib[i]:+.4f}  ({100*contrib[i]/primew:.1f}% of prime(w))")
    # is the defect μ from a few primes or spread?
    pos=contrib[contrib>0].sum(); neg=contrib[contrib<0].sum()
    print(f"  sum positive contribs={pos:.4f}, negative={neg:.4f}; primes 2,3 alone: {contrib[nn==2][0]:+.4f},{contrib[nn==3][0]:+.4f}")
print("="*90)
print("If a FEW low primes dominate prime(w) → the razor's arithmetic is controlled by the first primes (finite input).")
print("="*90)

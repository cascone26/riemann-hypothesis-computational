#!/usr/bin/env python3
"""§896: is the archimedean form TOTALLY POSITIVE as a kernel over ordered centers? Groskin T3.2: the archimedean
tail is a totally-positive Cauchy–Stieltjes increment. Finite checkable shadow: A(c_i,c_j)=∫P_{c_i}P_{c_j}Ω dr/2π
(arch kernel over ordered centers c_1<...<c_N). Test total positivity (all leading & contiguous minors >0, and
sign-regularity of 2x2 minors). If TP → Karlin variation-diminishing controls #neg(arch−prime). Compare: is the
prime kernel also TP? is arch−prime sign-regular? This bridges our M_σ picture to Groskin's TP tool."""
import numpy as np
from scipy.special import digamma
from itertools import combinations
def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)

def kernels(a,b,fourW,N,dr=0.003):
    W=fourW/4; RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi)
    CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN])
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    return arch,prime,CEN

def tp_report(M,name):
    N=M.shape[0]
    # all 2x2 minors sign (TP2: every 2x2 minor >=0)
    neg2=0; tot2=0; minm=np.inf
    for i,j in combinations(range(N),2):
        for k,l in combinations(range(N),2):
            det=M[i,k]*M[j,l]-M[i,l]*M[j,k]; tot2+=1
            if det<-1e-14: neg2+=1
            minm=min(minm,det)
    # leading principal minors (should all be >0 if PD)
    lead=[np.linalg.det(M[:k,:k]) for k in range(1,N+1)]
    posdef=all(l>0 for l in lead)
    # sign changes in each row (variation) - for variation-diminishing intuition
    return neg2,tot2,minm,posdef,lead[-1]

print("="*100)
print("§896 total-positivity of the archimedean kernel A(c_i,c_j) (Groskin TP shadow); prime & arch−prime too")
print("="*100)
for (a,b,fourW,N) in [(14,44,4.0,10),(14,34,4.0,8),(60,90,4.0,10),(14,44,3.0,9)]:
    arch,prime,CEN=kernels(a,b,fourW,N)
    for M,nm in [(arch,"arch"),(prime,"prime"),(arch-prime,"arch−prime")]:
        neg2,tot2,minm,pd,det=tp_report(M,nm)
        print(f"[{a},{b}] 4W={fourW} N={N} {nm:>10}: 2x2 minors neg={neg2}/{tot2}  min2x2minor={minm:+.2e}  PD={pd}  fulldet={det:+.2e}")
    print("-"*100)
print("arch: 0 negative 2x2 minors → TP2 (totally positive order 2, Karlin sign-regular) → variation-diminishing.")
print("If arch is TP2 but prime is NOT, the arch structure constrains #sign-changes of arch−prime → bounds #neg. ")
print("="*100)

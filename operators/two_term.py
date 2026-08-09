#!/usr/bin/env python3
"""§922: THE TWO-TERM ATTACK (survey-recommended). Prove arch+pole ⪰ p=2 term unconditionally. First verify it's
COMFORTABLE (not razor) across the low zone + bandwidths, and map the razor ONSET (which prime kills comfort).
margin_K = λmin(arch+pole − primes{n≤K}, G). K=2 (p=2 only) should be comfortable; K=3 razor onset."""
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
def build(a,b,fourW,frac=0.7,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    return P,Gm,arch,pole,RG,N,fourW
def pmat(P,RG,K,fourW,dr=0.004):
    N=P.shape[0]; m=(LOGN<fourW)&(NS<=K); ln=LOGN[m]; lw=LAMW[m]
    if len(ln)==0: return np.zeros((N,N))
    cu=np.cos(np.outer(ln,RG)); pr=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); pr[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    return 0.5*(pr+pr.T)
print("="*92)
print("§922 two-term margin: λmin(arch+pole − primes{n≤K}) across low-zone windows & bandwidths. K=2 comfortable?")
print("="*92)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'K=2(p2)':>10} {'K=3':>10} {'K=4':>10} {'K=5':>10} {'K=full':>11}")
for (a,b) in [(14,30),(14,44),(20,50),(30,60),(45,75)]:
    for fourW in [2.0,3.0,4.0]:
        P,Gm,arch,pole,RG,N,fW=build(a,b,fourW)
        row=[]
        for K in [2,3,4,5,3000]:
            pr=pmat(P,RG,K,fW); row.append(float(eigh(arch+pole-pr,Gm,eigvals_only=True)[0]))
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {row[0]:>+10.4f} {row[1]:>+10.4f} {row[2]:>+10.4f} {row[3]:>+10.4f} {row[4]:>+11.2e}")
print("="*92)
print("K=2 (arch+pole vs p=2 only): if uniformly comfortable (≫0) across the low zone → arb/SOS-certifiable unconditional")
print("result (arch dominates the first prime). Razor onset = the K where margin crashes to ~1e-3 or negative.")
print("="*92)

#!/usr/bin/env python3
"""§924: is pnorm(p=2)=λmax(p2 form,G) BOUNDED as bandwidth 4W→∞? (the load-bearing fact for proving arch+pole⪰p=2).
The p=2 term includes prime-POWERS of 2 (2,4,8,16,32,...) with log<4W. Compute pnorm2 vs 4W (large), and floor, to
find the explicit bound C2 and the height threshold where floor(t)>C2."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def build(a,b,fourW,frac=0.7,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    # p=2 powers: n=2,4,8,16,... with log n < 4W ; Lambda(2^k)=log2, weight log2/sqrt(2^k)
    p2=np.zeros((N,N)); k=1
    while k*np.log(2)<fourW:
        n=2**k; w=np.log(2.0)/np.sqrt(n); cu=np.cos(RG*np.log(n))
        for j in range(N):
            p2[j,:]+=2.0*w*((P[j][None,:]*P)@cu)*dr/(2*np.pi)
        k+=1
    p2=0.5*(p2+p2.T)
    floor=float(eigh(arch,Gm,eigvals_only=True)[0])  # arch floor (no pole, the pure archimedean min)
    pnorm2=float(eigh(p2,Gm,eigvals_only=True)[-1])
    return N,floor,pnorm2
print("="*80)
print("§924 is pnorm(p=2 powers) BOUNDED as 4W grows? floor(window) vs pnorm2. threshold floor>C2?")
print("="*80)
print(f"{'window':>12} {'4W':>5} {'N':>3} {'#2-powers':>9} {'floor(arch)':>11} {'pnorm(p2)':>10}")
for (a,b) in [(14,44),(20,50),(30,60)]:
    for fourW in [2.0,4.0,6.0,8.0,10.0,12.0]:
        N,floor,pnorm2=build(a,b,fourW)
        n2=int(fourW/np.log(2))
        print(f"{f'[{a},{b}]':>12} {fourW:>5.1f} {N:>3} {n2:>9} {floor:>11.4f} {pnorm2:>10.4f}")
    print("-"*80)
print("="*80)
print("If pnorm(p2) SATURATES to a constant C2 as 4W→∞ (bounded), and floor(t)=ln(t/2π)/2π grows past C2 → arch⪰p2")
print("for all t above an EXPLICIT threshold t2=2π·exp(2π·C2). The bounded C2 is the whole point (p=2 doesn't grow).")
print("="*80)

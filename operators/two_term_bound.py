#!/usr/bin/env python3
"""§923: PROVE arch+pole ⪰ p=2 via the operator bound margin ≥ floor − pnorm(p=2). Unlike the full prime form (pnorm
grows with bandwidth → §875 fails at corner), the SINGLE p=2 term is BOUNDED. Compute floor=λmin(arch+pole,G) and
pnorm2=λmax(p2form,G) across the low zone + bandwidths. If floor > pnorm2 UNIFORMLY → arch+pole ⪰ p=2 PROVEN (both
explicit/zero-free/arb-certifiable). Also report the true margin (arch+pole−p2) to see bound tightness."""
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
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    # p=2 term: 2*(log2/sqrt2)*(1/2pi) integral phi^2 cos(r log2)
    l2=np.log(2.0); w2=l2/np.sqrt(2.0)
    cu=np.cos(RG*l2); p2=np.zeros((N,N))
    for j in range(N):
        p2[j,:]=2.0*w2*((P[j][None,:]*P)@cu)*dr/(2*np.pi)
    p2=0.5*(p2+p2.T)
    floor=float(eigh(arch+pole,Gm,eigvals_only=True)[0])
    pnorm2=float(eigh(p2,Gm,eigvals_only=True)[-1])
    margin=float(eigh(arch+pole-p2,Gm,eigvals_only=True)[0])
    return N,floor,pnorm2,margin
print("="*96)
print("§923 PROVE arch+pole ⪰ p=2: is floor > pnorm(p=2) uniformly? (both explicit/zero-free → certifiable)")
print("="*96)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'floor':>9} {'pnorm(p2)':>10} {'floor−pnorm2':>12} {'true margin':>12} {'bound OK?':>9}")
allok=True; minbound=99
for (a,b) in [(14,30),(14,44),(20,50),(30,60),(45,75),(60,90)]:
    for fourW in [2.0,3.0,4.0,5.0]:
        N,floor,pnorm2,margin=build(a,b,fourW)
        bd=floor-pnorm2; ok=bd>0; allok=allok and ok; minbound=min(minbound,bd)
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {floor:>9.4f} {pnorm2:>10.4f} {bd:>+12.4f} {margin:>+12.4f} {'YES' if ok else 'NO':>9}")
print("="*96)
print(f"floor > pnorm(p=2) everywhere: {allok}   min(floor−pnorm2) = {minbound:+.4f}")
print("If YES uniformly → arch+pole ⪰ p=2 is PROVEN by the operator bound (margin ≥ floor−pnorm2 > 0), arb-certifiable.")
print("="*96)

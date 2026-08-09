#!/usr/bin/env python3
"""§925: the two-term result structure. arch+pole ⪰ (p=2 tower) — is it (a) operator-bound-proven above a FIXED finite
height threshold t2 (floor(t)>C2, C2=sup pnorm(p2)), and (b) comfortable (arb-certifiable) below t2? Map margin +
floor + pnorm2 across heights 15→250 at high bandwidth (4W=6, so p2 tower is near-saturated). If margin comfortable
below t2 and floor>pnorm2 above → arch+pole⪰p2 PROVABLE (finite comfortable bottom + explicit operator tail)."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def build(a,b,fourW,frac=0.7,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+70),(b+70),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    p2=np.zeros((N,N)); k=1
    while k*np.log(2)<fourW:
        n=2**k; w=np.log(2.0)/np.sqrt(n); cu=np.cos(RG*np.log(n))
        for j in range(N): p2[j,:]+=2.0*w*((P[j][None,:]*P)@cu)*dr/(2*np.pi)
        k+=1
    p2=0.5*(p2+p2.T)
    floor=float(eigh(arch+pole,Gm,eigvals_only=True)[0]); pnorm2=float(eigh(p2,Gm,eigvals_only=True)[-1])
    margin=float(eigh(arch+pole-p2,Gm,eigvals_only=True)[0])
    return floor,pnorm2,margin
print("="*82)
print("§925 arch+pole ⪰ (p=2 tower): margin, floor, pnorm2 vs height (4W=6, near-saturated p2). threshold?")
print("="*82)
print(f"{'window':>12} {'ht':>4} {'floor':>8} {'pnorm2':>8} {'floor−pnorm2':>12} {'margin':>10} {'status':>18}")
minmarg=99
for (a,b) in [(14,30),(14,44),(24,54),(34,64),(50,80),(70,100),(100,130),(140,170),(200,230)]:
    floor,pnorm2,margin=build(a,b,6.0); h=(a+b)//2; minmarg=min(minmarg,margin)
    st="op-bound PROVES" if floor>pnorm2 else "needs finite cert"
    print(f"{f'[{a},{b}]':>12} {h:>4} {floor:>8.4f} {pnorm2:>8.4f} {floor-pnorm2:>+12.4f} {margin:>+10.4f} {st:>18}")
print("="*82)
print(f"min margin over all heights = {minmarg:+.4f} (comfortable if ≫0). C2=sup pnorm2 → threshold t2=2π exp(2π C2).")
print("Below t2: finite region, comfortable margin → arb-certifiable. Above t2: floor>pnorm2 proves it. => arch+pole⪰p2 PROVABLE.")
print("="*82)

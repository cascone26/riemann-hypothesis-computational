#!/usr/bin/env python3
"""ATTACK the boundary-mode positivity crash (§813): §812 localized the uniform difficulty to the
boundary/corner-adjacent concentration modes. Decompose Q = pole+arch-prime on EACH concentration
eigenvector (core vs boundary) to find WHICH term crashes positivity on the boundary — arithmetic
(prime) or archimedean (arch). That pins the B18 target to a mechanism. 4W=1.46, 8-probe basis."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec=200
fourW=1.46; nprobes=8
ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W; Wa=ac.W
CEN=[round(14+30*i/(nprobes-1),3) for i in range(nprobes)]; n=len(CEN); cs=[arb(c) for c in CEN]
def sinc_sq(y):
    ay=abs(y); ymax=float(ay.mid())+float(ay.rad())
    if ymax<0.5:
        y2=y*y; return acb(1)-y2/3+2*y2*y2/45+acb(arb(0,ymax**6/45+1e-30),arb(0,ymax**6/45+1e-30))
    return (y.sin()/y)**2
def Fev(c,r): return sinc_sq(Wa*(r-acb(c)))+sinc_sq(Wa*(r+acb(c)))
splits=sorted(set([14.0,44.0]+[c for c in CEN if 14<c<44]))
# component matrices: pole, arch, prime, and G, A(concentration)
POLE=[[None]*n for _ in range(n)]; ARCH=[[None]*n for _ in range(n)]; PRIME=[[None]*n for _ in range(n)]
G=[[None]*n for _ in range(n)]; A=[[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j,n):
        POLE[j][k]=POLE[k][j]=ac.pole_entry(cs[j],cs[k])
        ARCH[j][k]=ARCH[k][j]=zc.arch_entry_v3(cs[j],cs[k],arb(400))
        PRIME[j][k]=PRIME[k][j]=ac.prime_entry(cs[j],cs[k])
        G[j][k]=G[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
        I=acb(0)
        for a,b in zip(splits[:-1],splits[1:]):
            m=max(2,int((b-a)))
            for tt in range(m): I=I+acb.integral(lambda r,_: Fev(CEN[j],r)*Fev(CEN[k],r), acb(a+(b-a)*tt/m),acb(a+(b-a)*(tt+1)/m))
        A[j][k]=A[k][j]=2*I.real
Gm=np.array([[float(G[j][k].mid()) for k in range(n)] for j in range(n)])
Am=np.array([[float(A[j][k].mid()) for k in range(n)] for j in range(n)])
Pm=np.array([[float(POLE[j][k].mid()) for k in range(n)] for j in range(n)])
ARm=np.array([[float(ARCH[j][k].mid()) for k in range(n)] for j in range(n)])
PRm=np.array([[float(PRIME[j][k].mid()) for k in range(n)] for j in range(n)])
gw,gv=np.linalg.eigh(Gm); Bf=gv/np.sqrt(gw)
ew,evc=np.linalg.eigh(Bf.T@Am@Bf); order=np.argsort(ew)[::-1]
print(f"4W={fourW}: decompose Q=pole+arch-prime on each concentration eigenvector (unit ‖v‖_G)")
print(" mode  conc     Q        pole     arch     prime    | in-window? crash?")
for i in range(n):
    v=Bf@evc[:,order[i]]                    # G-orthonormal eigenvector
    def quad(M): return float(v@np.array([[float(M[a][b].mid()) for b in range(n)] for a in range(n)])@v) if False else float(v@M@v)
    po=v@Pm@v; ar=v@ARm@v; pr=v@PRm@v; Q=po+ar-pr; nrm=v@Gm@v
    po,ar,pr,Q=po/nrm,ar/nrm,pr/nrm,Q/nrm
    tag="CORE" if ew[order[i]]>0.95 else ("BOUNDARY" if ew[order[i]]>0.2 else "out")
    print(f"  {i+1}   {ew[order[i]]:.3f}  {Q:+.4f}  {po:+.4f}  {ar:+.4f}  {pr:+.4f}  | {tag}")
print("\nREAD: compare CORE vs BOUNDARY rows — which term (arch or prime) collapses Q on the boundary?")

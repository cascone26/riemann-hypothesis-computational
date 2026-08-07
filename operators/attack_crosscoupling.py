#!/usr/bin/env python3
"""§813b: the §812 crash is CROSS-COUPLING (Q positive per-mode, §813). Find the low-eig direction of
the pencil (Q,G) on the top-6 concentration eigenspace, decompose it, and show WHICH mode-pair coupling
+ WHICH term (arch/prime) creates it. Reuses the component matrices at 4W=1.46."""
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
ARCH=np.zeros((n,n)); PRIME=np.zeros((n,n)); POLE=np.zeros((n,n)); G=np.zeros((n,n)); A=np.zeros((n,n))
for j in range(n):
    for k in range(j,n):
        POLE[j,k]=POLE[k,j]=float(ac.pole_entry(cs[j],cs[k]).mid())
        ARCH[j,k]=ARCH[k,j]=float(zc.arch_entry_v3(cs[j],cs[k],arb(400)).mid())
        PRIME[j,k]=PRIME[k,j]=float(ac.prime_entry(cs[j],cs[k]).mid())
        G[j,k]=G[k,j]=float((2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))).mid())
        I=acb(0)
        for a,b in zip(splits[:-1],splits[1:]):
            m=max(2,int((b-a)))
            for tt in range(m): I=I+acb.integral(lambda r,_: Fev(CEN[j],r)*Fev(CEN[k],r), acb(a+(b-a)*tt/m),acb(a+(b-a)*(tt+1)/m))
        A[j,k]=A[k,j]=2*float(I.real.mid())
Q=POLE+ARCH-PRIME
gw,gv=np.linalg.eigh(G); Bf=gv/np.sqrt(gw)
ew,evc=np.linalg.eigh(Bf.T@A@Bf); order=np.argsort(ew)[::-1]
K=6; V=Bf@evc[:,order[:K]]                       # top-6 concentration modes, G-orthonormal
# pencil (Q,G) on top-6: generalized eig
QE=V.T@Q@V; GE=V.T@G@V
w,U=np.linalg.eigh(np.linalg.solve(GE,QE))       # generalized eigenvalues/vecs
imin=np.argmin(w.real)
c=U[:,imin].real; c=c/np.sqrt(c@GE@c)            # the low-eig direction, in eigenmode coords
print(f"4W={fourW}: pencil λ_min on top-6 = {w.real.min():.5f}  (the §812 m_E)")
print(f"low-eig direction = combination of concentration modes (coeffs, modes 1..6):")
print("  " + "  ".join(f"m{i+1}:{c[i]:+.3f}" for i in range(K)))
# which modes dominate the low-eig direction?
dom=np.argsort(-np.abs(c))[:3]
print(f"  dominant modes: {[f'm{d+1}(conc{ew[order[d]]:.2f})' for d in dom]}")
# decompose Q on this direction into arch/prime/pole
vphys=V@c
po=vphys@POLE@vphys; ar=vphys@ARCH@vphys; pr=vphys@PRIME@vphys; nrm=vphys@G@vphys
print(f"\nQ on the low-eig direction: pole={po/nrm:+.4f} arch={ar/nrm:+.4f} prime={pr/nrm:+.4f} → Q={(po+ar-pr)/nrm:+.5f}")
print("compare to per-mode arch~+0.23: is arch SUPPRESSED or prime AMPLIFIED on this direction?")
# also: arch and prime pencil min-eigs separately (which term drives the low-eig?)
for name,M in [("arch",ARCH),("prime",PRIME),("Q",Q)]:
    ME=V.T@M@V; wm=np.linalg.eigvalsh(np.linalg.solve(GE,ME) if False else np.linalg.inv(GE)@ME)
    print(f"  {name}: pencil eigs on top-6 = {np.sort(wm.real).round(3)}")

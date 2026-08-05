import sys, math
sys.path.insert(0,'.')
from arb_ldlt_certify import build_arb_tau
from flint import arb, ctx
import mpmath as mp, numpy as np
def spectrum(c,N,prec):
    ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec)
    mp.mp.dps=int(prec*0.29)
    Amp=mp.matrix(DIM,DIM)
    for i in range(DIM):
        for j in range(DIM): Amp[i,j]=mp.mpf(A[i,j].str(int(prec*0.28),radius=False))
    E=sorted(float(mp.log10(abs(e))) for e in mp.eigsy(Amp,eigvals_only=True))
    return E  # log10 of |eigenvalues|, ascending
print("A) true MIN-EIG vs N (c=30):  is log10(min-eig) linear or quadratic in N?")
Ns=[]; mins=[]
for N in [4,6,8,10,12]:
    E=spectrum(30,N,1400); Ns.append(N); mins.append(E[0])
    print(f"  N={N:2} dim={2*N+1:2} log10(min-eig)={E[0]:+.2f}  cascade(lowest5)={[f'{x:+.1f}' for x in E[:5]]}")
lin=np.polyfit(Ns,mins,1); quad=np.polyfit(Ns,mins,2)
print(f"  linear: {lin[0]:.3f}N{lin[1]:+.2f} (resid {np.std(np.polyval(lin,Ns)-mins):.2f}) | quad: {quad[0]:+.4f}N^2{quad[1]:+.3f}N{quad[2]:+.1f} (resid {np.std(np.polyval(quad,Ns)-mins):.3f})")
print("B) true MIN-EIG vs L=log c (N=8):")
Ls=[]; m2=[]
for c in [13,30,100,300,1000]:
    E=spectrum(c,8,1400); L=math.log(c); Ls.append(L); m2.append(E[0])
    print(f"  c={c:5} L={L:.3f} log10(min-eig)={E[0]:+.2f}")
lin2=np.polyfit(Ls,m2,1)
print(f"  linear in L: {lin2[0]:.3f}*L{lin2[1]:+.2f} (resid {np.std(np.polyval(lin2,Ls)-m2):.2f})")
print("C) cascade geometry (c=30,N=10): consecutive gaps in log10(eig):")
E=spectrum(30,10,1600); gaps=[E[i+1]-E[i] for i in range(len(E)-1)]
print("  ", [f'{g:.1f}' for g in gaps[:12]])

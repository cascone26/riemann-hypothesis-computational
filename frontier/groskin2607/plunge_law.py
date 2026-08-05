import sys, math, numpy as np
sys.path.insert(0,'.')
from arb_ldlt_certify import build_arb_tau
from flint import arb, ctx
import mpmath as mp
def mineig(c,N):
    prec=max(700,int((2.6*N+120)*4.2)); ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec); mp.mp.dps=int(prec*0.29)
    Amp=mp.matrix(DIM,DIM)
    for i in range(DIM):
        for j in range(DIM): Amp[i,j]=mp.mpf(A[i,j].str(int(prec*0.28),radius=False))
    return min(float(mp.log10(abs(e))) for e in mp.eigsy(Amp,eigvals_only=True))
# 1) pin the N-law exponent: fit log10|mineig| = a + b*N^s, find best s (c=30)
Ns=np.array([8,12,16,20,24,28,32,36,40,44]); y=np.array([-26.73,-35.52,-43.72,-50.31,-56.78,-63.15,-68.36,-73.26,-78.03,-82.50])
best=None
for s in np.linspace(0.3,1.0,71):
    X=np.vstack([np.ones_like(Ns,dtype=float),Ns.astype(float)**s]).T
    cc,_,_,_=np.linalg.lstsq(X,y,rcond=None); r=np.std(X@cc-y)
    if best is None or r<best[0]: best=(r,s,cc[1])
print(f"N-law (c=30): best log10|mineig| ~ {best[2]:.2f}*N^{best[1]:.2f}  (resid {best[0]:.3f})   [s=0.5 => prolate-type sqrt-N plunge]")
# 2) min-eig vs c at fixed N=28 (does it decrease steadily with c = band growing?)
print("\nmin-eig vs c at fixed N=28 (band resolution fixed):")
res=[]
for c in [5,7,13,30,70,150,400,1000]:
    v=mineig(c,28); L=math.log(c); res.append((c,L,v))
    print(f"  c={c:5} L=log c={L:.3f}  log10(min-eig)={v:+.2f}")
Ls=np.array([r[1] for r in res]); ys=np.array([r[2] for r in res])
m=np.polyfit(Ls,ys,1); print(f"  fit: log10(min-eig) ~ {m[0]:.2f}*L {m[1]:+.2f}  => min-eig ~ c^{m[0]*2.303:.1f}")

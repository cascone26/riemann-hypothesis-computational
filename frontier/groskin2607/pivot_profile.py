import sys, math
sys.path.insert(0,'.')
from arb_ldlt_certify import build_arb_tau
from flint import arb, ctx
import numpy as np
def pivots_and_eig(c,N,prec):
    ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec)
    # pivots
    d=[None]*DIM; Lf=[[arb(0)]*DIM for _ in range(DIM)]
    for i in range(DIM):
        s=A[i,i]
        for k in range(i): s=s-Lf[i][k]*Lf[i][k]*d[k]
        d[i]=s
        for j in range(i+1,DIM):
            t=A[j,i]
            for k in range(i): t=t-Lf[j][k]*Lf[i][k]*d[k]
            Lf[j][i]=t/d[i]
    piv=[float(x.mid() if hasattr(x,'mid') else x) for x in d]
    # high-prec float eigenvalues via mpmath (proper margin)
    import mpmath as mp
    mp.mp.dps=int(prec*0.3)
    Amp=mp.matrix(DIM,DIM)
    for i in range(DIM):
        for j in range(DIM):
            Amp[i,j]=mp.mpf(A[i,j].str(int(prec*0.29), radius=False))
    E=mp.eigsy(Amp, eigvals_only=True)
    eigs=sorted(float(e) for e in E)
    return piv,eigs
c,N=30,8
piv,eigs=pivots_and_eig(c,N,500)
print(f"c={c} N={N} dim={2*N+1}")
print("pivots by freq index n=-N..N:")
for i,p in enumerate(piv): print(f"  n={i-N:+3d}: pivot={p:+.4e}  log10|{math.log10(abs(p)):+.2f}|")
print(f"min pivot = {min(piv,key=abs):.4e} at n={piv.index(min(piv,key=abs))-N}")
print("\ntrue eigenvalues (mpmath eigsy, sorted):")
for e in eigs[:5]: print(f"  {e:+.5e}  log10={math.log10(abs(e)):+.2f}")
print(f"  ... max={eigs[-1]:.4e}")
print(f"\nMIN EIGENVALUE = {eigs[0]:.4e} (log10={math.log10(abs(eigs[0])):+.2f})  vs MIN PIVOT log10={math.log10(abs(min(piv,key=abs))):+.2f}")
print(f"count eigenvalues below 1e-10: {sum(1 for e in eigs if abs(e)<1e-10)} / {len(eigs)}")

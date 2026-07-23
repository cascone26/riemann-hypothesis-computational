#!/usr/bin/env python3
"""§198: extract the explicit 4×4 even-mode Weil matrix (modes k=1,3,5,7) + its arch/prime split + Sylvester minors,
at a grid across (0.42, ½log3). Goal: see the entries' a-dependence (recognize closed forms) and confirm all 4
leading principal minors of the FULL even-block stay >0 (⟺ PD ⟺ positivity). mpmath; Arb-certify next."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices
mp.mp.dps=40; K=24; M=800
ev=[c for c in range(K) if (c+1)%2==1][:4]   # k=1,3,5,7
def minors(Qf,Mass):
    # generalized: minors of L^{-1}Qf L^{-T} where Mass=LL^T (Cholesky) -> use Mass-normalized block
    n=4; A=mp.matrix(n,n); Me=mp.matrix(n,n)
    for i in range(n):
        for j in range(n): A[i,j]=Qf[ev[i],ev[j]]; Me[i,j]=Mass[ev[i],ev[j]]
    # whiten by Mass: B = Me^{-1/2} A Me^{-1/2}; PD(A wrt Me) <=> PD(B). Use eig of (A,Me) leading minors via chol.
    Lc=mp.cholesky(Me); Linv=Lc**-1
    B=Linv*A*Linv.T
    ms=[]
    for k in range(1,n+1):
        sub=mp.matrix(k,k)
        for i in range(k):
            for j in range(k): sub[i,j]=B[i,j]
        ms.append(mp.det(sub))
    return ms, A, Me
print(f"§198 4×4 even Weil block (k=1,3,5,7)  M={M}")
for a in [0.44,0.48,0.52,0.549]:
    am=mp.mpf(str(a)); Qf,Qa,Mass,Vv,xs,w=weil_matrices(am,K,M)
    ms,A,Me=minors(Qf,Mass)
    allpos=all(m>0 for m in ms)
    print(f"a={a:.3f} Sylvester minors: {[mp.nstr(m,4) for m in ms]}  {'ALL>0 (PD✓)' if allpos else 'FAIL'}")
print("\nREAD: all 4 minors >0 at every a ⟹ the 4×4 is PD ⟹ positivity holds (mod tail). The minors' a-shapes guide the")
print("closed-form/Arb-certified proof. det(4×4)→0-ish tracks the near-null (last minor is the tightest).")

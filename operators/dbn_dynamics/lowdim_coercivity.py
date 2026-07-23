#!/usr/bin/env python3
"""§197: is the even-sector coercivity effectively LOW-dimensional? Near-null ≈ 0.96·cos + (−0.28)·mode3 + tiny
(§195) ⟹ maybe 2-3 odd-k modes capture the margin. Compute Weil-form λ_min restricted to the first d odd-k modes
(d=1,2,3,4) vs the full even margin. If d=2 or 3 reproduces it, the crux is a SMALL explicit matrix positivity."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=800
ev=[c for c in range(K) if (c+1)%2==1]   # odd-k modes: k=1,3,5,...
def block(Q,idx):
    B=mp.matrix(len(idx),len(idx))
    for i in range(len(idx)):
        for j in range(len(idx)): B[i,j]=Q[idx[i],idx[j]]
    return B
print(f"§197 low-dim even coercivity  K={K} M={M}  (odd-k modes k=1,3,5,7,...)")
print(f"{'a':>7} {'λmin[d=1]':>12} {'d=2':>12} {'d=3':>12} {'d=4':>12} {'full(d=12)':>12}")
for a in [0.44,0.48,0.52,0.549]:
    am=mp.mpf(str(a)); Qf,Qa,Mass,Vv,xs,w=weil_matrices(am,K,M)
    row=[]
    for d in [1,2,3,4,len(ev)]:
        idx=ev[:d]
        E,_=gen_eig_sym(block(Qf,idx),block(Mass,idx))
        row.append(mp.nstr(E[0],5))
    print(f"{a:>7.3f} {row[0]:>12} {row[1]:>12} {row[2]:>12} {row[3]:>12} {row[4]:>12}",flush=True)
print("\nREAD: the smallest d whose λmin ≈ full margin = the true dimension of the crux. d=2/3 ⟹ a 2x2/3x3 explicit")
print("matrix positivity (digamma entries) — a FINITE, hand-attackable inequality. d large ⟹ genuinely high-dim coercivity.")

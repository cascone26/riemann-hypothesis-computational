#!/usr/bin/env python3
"""§195: is the even near-null a FIXED function in ξ=x/a? Basis sin(kπ(ξ+1)/2) is a-independent, so fixed shape
⟺ the even ground-state COEFFICIENT vector (Euclidean-normalized) is constant across a. Compute cross-a dot
products. If ≈1, margin(a) is a FIXED function's Rayleigh quotient under the a-varying Weil operator."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=600
evidx=[c for c in range(K) if (c+1)%2==1]
def sub(Q,idx):
    n=len(idx); B=mp.matrix(n,n)
    for i in range(n):
        for j in range(n): B[i,j]=Q[idx[i],idx[j]]
    return B
vecs=[]; alist=[0.44,0.48,0.52,0.549]
for a in alist:
    am=mp.mpf(str(a)); Qf,Qa,Mass,Vv,xs,w=weil_matrices(am,K,M)
    Ef,Vf=gen_eig_sym(sub(Qf,evidx),sub(Mass,evidx))
    c=[Vf[r,0] for r in range(len(evidx))]
    nrm=mp.sqrt(mp.fsum(x*x for x in c)); c=[x/nrm for x in c]
    if c[0]<0: c=[-x for x in c]                       # sign-fix
    vecs.append(c)
print(f"§195 fixed-shape (even near-null coeff vector, Euclidean-normalized) K={K}")
print("cross-a |dot| (1.0 = identical shape in ξ=x/a):")
for i in range(len(alist)):
    row=[mp.nstr(abs(mp.fsum(vecs[i][t]*vecs[j][t] for t in range(len(evidx)))),6) for j in range(len(alist))]
    print(f"  a={alist[i]:.3f}: "+"  ".join(row))
print("\ntop-5 coeffs (odd-k sines) at each a — are they the SAME numbers?")
for i,a in enumerate(alist):
    print(f"  a={a:.3f}: "+"  ".join(mp.nstr(vecs[i][t],4) for t in range(5)))
print("\nREAD: dot≈1 across a ⟹ FIXED near-null function in ξ; then margin(a)=⟨v_fixed,(Qa−P2)(a)v_fixed⟩ — a SINGLE")
print("explicit function's a-varying Weil form, the cleanest possible proof object.")

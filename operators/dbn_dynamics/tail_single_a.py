#!/usr/bin/env python3
"""Test tail convergence for a single critical a value."""
import sys
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=800

ev=[c for c in range(K) if (c+1)%2==1]

def block(Q,idx):
    B=mp.matrix(len(idx),len(idx))
    for i in range(len(idx)):
        for j in range(len(idx)): B[i,j]=Q[idx[i],idx[j]]
    return B

a_test=mp.mpf('0.549')
print(f"Testing a={a_test}...")
Qf,Qa,Mass,Vv,xs,w=weil_matrices(a_test,K,M)

print("\nDimension convergence for a=0.549:")
print(f"{'d':>3} {'λmin(d)':>18} {'Tail from full':>18}")
print("-"*45)

for d in [1,2,3,4,6,8,10,12]:
    if d <= len(ev):
        idx=ev[:d]
        E,_=gen_eig_sym(block(Qf,idx),block(Mass,idx))
        lam_d=E[0]

        if d==len(ev):
            lam_full=lam_d
            tail=mp.mpf(0)
        else:
            tail=None

        if tail is not None:
            print(f"{d:>3} {mp.nstr(lam_d,15):>18} {mp.nstr(tail,15):>18}")
        else:
            print(f"{d:>3} {mp.nstr(lam_d,15):>18}  (full value)")

# Now compute tail bounds
print("\n" + "="*45)
print("Computing full matrix for tail reference...")
idx_full=ev
E_full,_=gen_eig_sym(block(Qf,idx_full),block(Mass,idx_full))
lam_full=E_full[0]

print(f"\nFull λmin = {mp.nstr(lam_full,15)}")
print(f"\n{'d':>3} {'λmin(d)':>18} {'Tail(d)':>18} {'Ratio(tail/full)':>18}")
print("-"*60)
for d in [1,2,3,4,6,8,10,12]:
    if d <= len(ev):
        idx=ev[:d]
        E,_=gen_eig_sym(block(Qf,idx),block(Mass,idx))
        lam_d=E[0]
        tail=lam_d-lam_full
        if lam_full > 0:
            ratio=tail/lam_full
        else:
            ratio=mp.mpf('inf')
        print(f"{d:>3} {mp.nstr(lam_d,15):>18} {mp.nstr(tail,15):>18} {mp.nstr(ratio,15):>18}")

#!/usr/bin/env python3
"""Minimal version: just a=0.48 and a=0.549 with essential d values."""
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

print("§199 ESSENTIAL TAIL DATA (fast)")
print("="*90)

for a in [0.48, 0.549]:
    a_val=mp.mpf(str(a))
    Qf,Qa,Mass,Vv,xs,w=weil_matrices(a_val,K,M)

    # Get full
    idx_full=ev
    E_full,_=gen_eig_sym(block(Qf,idx_full),block(Mass,idx_full))
    lam_full=E_full[0]

    print(f"\n a = {a}:")
    print(f"  λmin(full) = {mp.nstr(lam_full,14)}")
    print()
    print(f"  {'d':>3}  {'λmin(d)':>15}  {'Tail':>15}  {'%':>10}")
    print("  " + "-"*58)

    for d in [4,6,8,10,12]:
        if d <= len(ev):
            idx=ev[:d]
            E,_=gen_eig_sym(block(Qf,idx),block(Mass,idx))
            lam_d=E[0]
            tail=lam_d-lam_full
            pct=(tail/lam_full)*100 if lam_full>0 else mp.mpf('inf')
            print(f"  {d:>3}  {mp.nstr(lam_d,12):>15}  {mp.nstr(tail,12):>15}  {mp.nstr(pct,8):>10}%")

#!/usr/bin/env python3
"""Fast tail convergence analysis: focus on critical test points only.
Uses mpmath but tests only the essential d values: 4,6,8,10,12."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=800

ev=[c for c in range(K) if (c+1)%2==1]   # odd-k modes

def block(Q,idx):
    B=mp.matrix(len(idx),len(idx))
    for i in range(len(idx)):
        for j in range(len(idx)): B[i,j]=Q[idx[i],idx[j]]
    return B

print("§199 FAST TAIL CONVERGENCE (key test points)")
print("="*80)
print(f"{'a':>7} {'d=4':>12} {'d=6':>12} {'d=8':>12} {'d=10':>12} {'d=12':>12} {'Tail(8)':>12} {'Tail(10)':>12}")
print("-"*100)

# Only the critical a values
for a in [0.45, 0.48, 0.50, 0.549]:
    am=mp.mpf(str(a))
    Qf,Qa,Mass,Vv,xs,w=weil_matrices(am,K,M)

    row=[mp.mpf(str(a))]
    results={}

    for d in [4,6,8,10,12]:
        if d <= len(ev):
            idx=ev[:d]
            E,_=gen_eig_sym(block(Qf,idx),block(Mass,idx))
            results[d]=E[0]

    idx_full=ev
    E_full,_=gen_eig_sym(block(Qf,idx_full),block(Mass,idx_full))
    results['full']=E_full[0]

    # Print eigenvalues for d=4,6,8,10,12
    for d in [4,6,8,10,12]:
        print(f"{mp.nstr(results[d],6):>12} ", end='')

    # Print tails (diff from full)
    if 8 in results:
        tail_8=results[8]-results['full']
        print(f"{mp.nstr(tail_8,6):>12} ", end='')
    else:
        print("      -      ", end='')

    if 10 in results:
        tail_10=results[10]-results['full']
        print(f"{mp.nstr(tail_10,6):>12}", end='')
    else:
        print("      -     ", end='')

    print()

print("-"*100)
print("\nINTERPRETATION:")
print("  Tail(8) = λmin(d=8) - λmin(full) ≈ impact of modes d>8")
print("  Tail(10) = λmin(d=10) - λmin(full) ≈ impact of modes d>10")
print("  If Tail(8) << full margin (< 1%), d=8 is safe for rigorous proof + tail bound.")

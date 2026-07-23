#!/usr/bin/env python3
"""§193: the analytic target's shape. Trace the CRITICAL EVEN sector across the p=2 window a∈(½log2,½log3):
α_even(a) (the single even arch negative eigenvalue) and full λmin_even(a) (the margin to prove ≥0). If these
have a clean/monotone form, the secular bound is a tractable 1-variable inequality. EVEN-v = odd-k sines."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=600
def evenblock(Q,K):
    idx=[c for c in range(K) if (c+1)%2==1]; n=len(idx); B=mp.matrix(n,n)
    for i in range(n):
        for j in range(n): B[i,j]=Q[idx[i],idx[j]]
    return B
l2,l3=0.5*math.log(2),0.5*math.log(3)
print(f"§193 even-sector critical margin  a∈(½log2={l2:.4f}, ½log3={l3:.4f})  K={K}")
print(f"{'a':>7} {'α_even (arch neg)':>18} {'margin λmin_even':>18} {'margin/|α_even|':>15}")
for a in [0.36,0.40,0.44,0.48,0.52,0.549]:
    am=mp.mpf(str(a)); Qf,Qa,Mass,Vv,xs,w=weil_matrices(am,K,M)
    Ea,_=gen_eig_sym(evenblock(Qa,K),evenblock(Mass,K))
    Ef,_=gen_eig_sym(evenblock(Qf,K),evenblock(Mass,K))
    r=Ef[0]/abs(Ea[0]) if Ea[0]<0 else mp.mpf('nan')
    print(f"{a:>7.3f} {mp.nstr(Ea[0],7):>18} {mp.nstr(Ef[0],7):>18} {mp.nstr(r,5):>15}",flush=True)
print("\nREAD: α_even(a) and margin(a) shapes. Clean/monotone ⟹ the secular bound is a tractable 1-var inequality;")
print("margin>0 across the whole window = the target. margin/|α_even| = the fractional lift the involution must supply.")

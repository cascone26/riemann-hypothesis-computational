#!/usr/bin/env python3
"""§194: WHY does the even-sector margin lock flat at ~2.6e-5 while α_even deepens 20×? Test: is the even ground
state a FIXED shape across a (⟹ fixed Rayleigh quotient)? Decompose the ground state v0 into its arch-Rayleigh
⟨v0,Qa v0⟩ and prime-Rayleigh ⟨v0,(Qf−Qa)v0⟩ parts, and measure the cross-a overlap of v0 (fixed shape?)."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=600
def evidx(K): return [c for c in range(K) if (c+1)%2==1]
def sub(Q,idx):
    n=len(idx); B=mp.matrix(n,n)
    for i in range(n):
        for j in range(n): B[i,j]=Q[idx[i],idx[j]]
    return B
prev=None
print(f"§194 flat-lock probe (even sector)  K={K}")
print(f"{'a':>7} {'margin':>12} {'arch-Rayleigh':>14} {'prime-Rayleigh':>15} {'overlap w/prev v0':>18}")
for a in [0.48,0.52,0.549]:
    am=mp.mpf(str(a)); Qf,Qa,Mass,Vv,xs,w=weil_matrices(am,K,M)
    idx=evidx(K); Qae=sub(Qa,idx); Qfe=sub(Qf,idx); Me=sub(Mass,idx)
    Ef,Vf=gen_eig_sym(Qfe,Me)
    v0=[Vf[r,0] for r in range(len(idx))]                 # ground state (Mass-normalized)
    def ray(Q): return mp.fsum(v0[i]*Q[i,j]*v0[j] for i in range(len(idx)) for j in range(len(idx)))
    aR=ray(Qae); pR=ray(Qfe)-aR
    ov='—'
    if prev is not None:
        ov=mp.nstr(abs(mp.fsum(v0[i]*Me[i,j]*prev[j] for i in range(len(idx)) for j in range(len(idx)))),6)
    print(f"{a:>7.3f} {mp.nstr(Ef[0],6):>12} {mp.nstr(aR,6):>14} {mp.nstr(pR,7):>15} {ov:>18}",flush=True)
    prev=v0
print("\nREAD: margin = arch-Rayleigh + prime-Rayleigh. If prime-Rayleigh ≈ −arch-Rayleigh + 2.6e-5 (they co-track) AND")
print("cross-a overlap≈1 (fixed ground-state shape), the flat-lock = a FIXED near-null vector whose arch/prime parts")
print("cancel to a constant — the finite-window analog of the CvS k_λ. That fixed vector is the proof handle.")

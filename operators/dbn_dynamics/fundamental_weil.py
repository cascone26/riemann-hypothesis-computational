#!/usr/bin/env python3
"""§196: Weil form of the EXACT fundamental even mode cos(πx/2a) = the k=1 basis function ⟹ just the (1,1) matrix
entry. Extract arch(cos)=Qa[0,0]/Mass[0,0], prime(cos)=(Qf−Qa)[0,0]/Mass[0,0], Weil(cos)=Qf[0,0]/Mass[0,0] across a.
Compare Weil(cos) to the true margin (§193): if close, the fundamental mode captures the crux; and arch(cos)(a),
prime(cos)(a) are now EXPLICIT closed-form-candidate functions of a to fit/prove."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=800
ev=[c for c in range(K) if (c+1)%2==1]; i0=ev[0]   # k=1 index
print(f"§196 fundamental-mode Weil form  cos(πx/2a)=basis k=1   K={K} M={M}")
print(f"{'a':>7} {'arch(cos)':>13} {'prime(cos)':>13} {'Weil(cos)':>13} {'true margin':>13} {'Weil(cos)/margin':>16}")
def evblk(Q):
    B=mp.matrix(len(ev),len(ev))
    for i in range(len(ev)):
        for j in range(len(ev)): B[i,j]=Q[ev[i],ev[j]]
    return B
for a in [0.44,0.48,0.52,0.549]:
    am=mp.mpf(str(a)); Qf,Qa,Mass,Vv,xs,w=weil_matrices(am,K,M)
    mass00=Mass[i0,i0]
    archc=Qa[i0,i0]/mass00; weilc=Qf[i0,i0]/mass00; primec=weilc-archc
    Ef,_=gen_eig_sym(evblk(Qf),evblk(Mass)); margin=Ef[0]
    print(f"{a:>7.3f} {mp.nstr(archc,6):>13} {mp.nstr(primec,6):>13} {mp.nstr(weilc,6):>13} {mp.nstr(margin,6):>13} {mp.nstr(weilc/margin,5):>16}",flush=True)
print("\nREAD: Weil(cos)=arch(cos)−|prime(cos)|. If Weil(cos) small-positive ≈ margin, the fundamental mode IS the near-")
print("minimizer & the crux is arch(cos)(a) > prime(cos)(a) — two explicit 1-var functions. If Weil(cos)≫margin, the")
print("drift correction (§195) matters and the true minimizer dips below the cosine — record which.")

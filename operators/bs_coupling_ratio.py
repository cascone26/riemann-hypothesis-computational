#!/usr/bin/env python3
"""Cycle 9 (2026-08-26): Birman-Schwinger coupling ratio rho(a) in the Suzuki screw frame.
Q = Qtail(+Qpole) - P,  A:=Qtail(+Qpole) arch-kinetic (PD),  P:=-(Qcore+Qprime) potential.
lambda_a>0  <=>  A - P > 0  <=>  rho(a):=sup_v P(v)/A(v) = lambda_max( A^{-1} P ) < 1.
KEY: rho(a) is an O(1) RATIO of two O(a) forms => numerically RESOLVABLE where lambda_a itself
(a difference of cancelling O(a) terms) sits below the float floor (cycle 8). If rho(a)<1 with a
uniform margin as a grows => Birman-Schwinger (brick 1758) gives lambda_a>0 forall a => RH.
Verify-before-ship: check A is PD (all eigs>0, well-conditioned) before trusting the ratio."""
import numpy as np, suzuki_exact_Q as EQ
from scipy.linalg import eigh
K=80
print(f"{'a':>5} {'rho(a)=supP/A':>13} {'margin 1-rho':>12} {'A mineig':>10} {'A cond':>10} {'#neg(A-P)':>10}")
rows=[]
for a in [0.5,0.8,1.0,1.3,1.6,2.0,2.5,3.0,3.5,4.0]:
    Q,pc=EQ.build_Q(a,K,return_pieces=True)
    def sym(M): return (M+M.T)/2
    A=sym(pc["tail"]+pc["pole"])       # arch kinetic (should be PD)
    P=sym(-(pc["core"]+pc["prime"]))   # potential (core is PSI*a*I<0 => -core>0; prime mostly<0)
    ea=eigh(A,eigvals_only=True)
    Amin=ea[0]; cond=ea[-1]/max(ea[0],1e-30)
    # generalized eig P v = rho A v  => rho_max
    try:
        rho=eigh(P,A,eigvals_only=True)[-1]
    except Exception as e:
        rho=float('nan')
    nneg=int((eigh(sym(Q),eigvals_only=True)<-1e-9).sum())
    rows.append((a,rho,Amin,cond))
    print(f"{a:5.1f} {rho:13.6f} {1-rho:12.6f} {Amin:10.3e} {cond:10.2e} {nneg:10d}")
print("\nrho<1 with margin & A PD (Amin>0) => coupling sub-unit => lambda_a>0 (brick 1758). Track margin vs a.")
# trend of margin
import numpy as np
aa=np.array([r[0] for r in rows]); mg=np.array([1-r[1] for r in rows])
print(f"margin(a): min={mg.min():.4f} at a={aa[np.argmin(mg)]:.1f}; trend last-3 slope={np.polyfit(aa[-3:],mg[-3:],1)[0]:+.4f}/a")

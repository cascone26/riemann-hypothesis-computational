import sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
sys.path.insert(0,'/Users/scones/projects/riemann/external/groskin_2607.02828/code')
from flint import ctx
import numpy as np
from arb_ldlt_certify import build_arb_tau
from dh_certified_tau import build_tau_general, dh_prime_side
from dh_exact import kj_exact
def spec(A,DIM):
    M=np.array([[float(A[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
    return np.linalg.eigvalsh(0.5*(M+M.T))
def mu2(ev):
    pos=ev[ev>1e-12]; d=len(pos); m1=pos.sum()/d; m2=(pos**2).sum()/d
    return m2/m1**2, d
print(f"{'c':>4} {'fn':>5} {'n_neg':>5} {'d_pos':>5} {'mu2':>8} {'2-mu2':>8}")
for c in [50,100]:
    N=64;prec=800;ctx.prec=prec
    Az,Dz=build_arb_tau(c,N,prec); ez=spec(Az,Dz)
    kj=float(kj_exact('0.75',5,c,prec).mid()); dw,dpos=dh_prime_side(c,prec)
    Ad,Dd=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=dw,prime_positions=dpos)
    ed=spec(Ad,Dd)
    for tag,ev in [('zeta',ez),('D-H',ed)]:
        m,d=mu2(ev); nn=int((ev<-1e-9).sum())
        print(f"{c:>4} {tag:>5} {nn:>5} {d:>5} {m:>8.4f} {2-m:>8.4f}",flush=True)
print("mu2 is SCALE-INVARIANT (mean-normalized 2nd moment). AF target for zeta ~ 4/3=1.333, proportion 2-mu2~0.667.")

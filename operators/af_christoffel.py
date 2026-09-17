import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import ctx
import numpy as np
from dh_certified_tau import build_arb_tau
def christoffel_bound(mu, m):
    # mu[0..2m] normalized moments (mu[0]=1). Lambda_m(0)=1/[H^{-1}]_00, bound = 1 - Lambda_m(0)
    H=np.array([[mu[i+j] for j in range(m+1)] for i in range(m+1)])
    Hinv=np.linalg.inv(H)
    Lam=1.0/Hinv[0,0]
    return 1.0-Lam
print("Christoffel-Markov proportion bound n+/d >= 1 - Lambda_m(0) from our EXACT certified moments:")
for c in [50,100]:
    ctx.prec=800
    A,DIM=build_arb_tau(c,64,800)
    M=np.array([[float(A[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
    ev=np.linalg.eigvalsh(0.5*(M+M.T)); ev=ev[ev>1e-12]; d=len(ev)
    mean=ev.sum()/d
    mu=[(ev**k).sum()/d/mean**k for k in range(5)]  # mean-normalized moments mu[0..4], mu[0]=1
    b1=christoffel_bound(mu,1); b2=christoffel_bound(mu,2)
    print(f"c={c}: mu=[1,{mu[1]:.3f},{mu[2]:.3f},{mu[3]:.3f},{mu[4]:.3f}]")
    print(f"      m=1 bound (moments<=k2): 1-Lam_1(0) = {b1:.4f}   (AF k=1 gives 2/3=0.667)")
    print(f"      m=2 bound (moments<=k4): 1-Lam_2(0) = {b2:.4f}   (AF HL*(4) target 13/18=0.722)")
# also target-moment bounds (RMT/HL exact): mu=1,1,4/3,2,13/4
mt=[1,1,4/3,2,13/4]
print(f"\nWith AF's TARGET moments mu=[1,1,4/3,2,13/4]: m=1 -> {christoffel_bound(mt,1):.4f}, m=2 -> {christoffel_bound(mt,2):.4f} (expect ~2/3, ~13/18)")

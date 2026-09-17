import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import ctx
import numpy as np
from dh_certified_tau import build_arb_tau
c=100
print(f"mu4 convergence toward HL target 13/4=3.25 (and mu2->4/3=1.333), c={c}, growing N:")
for N in [64,96,128]:
    ctx.prec=800
    A,DIM=build_arb_tau(c,N,800)
    M=np.array([[float(A[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
    ev=np.linalg.eigvalsh(0.5*(M+M.T)); ev=ev[ev>1e-12]; d=len(ev); mean=ev.sum()/d
    mu2=(ev**2).sum()/d/mean**2; mu4=(ev**4).sum()/d/mean**4
    print(f"  N={N:3d} dim={DIM} d(pos)={d}: mu2={mu2:.4f} (t 1.3333)  mu4={mu4:.4f} (t 3.2500)",flush=True)
print("mu4 climbing toward 3.25 with N would numerically realize AF's HL*(4) -> 13/18 on our certified instrument.")

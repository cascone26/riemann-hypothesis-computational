import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
import dh_certified_tau as DH
from dh_certified_tau import build_tau_general, margin, dh_prime_side
c,N,prec=14,48,800
dw,dpos=dh_prime_side(c,prec)
print(f"Locate the n_neg transition precisely (c={c},N={N}):")
for kj in [0.5940,0.5948,0.5950,0.5952,0.5954,0.5956,0.5958]:
    A,DIM=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=dw,prime_positions=dpos)
    nn,pd,mg=margin(A,DIM)
    print(f"  (k+J)={kj:.4f}: n_neg={nn}",flush=True)

import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import ctx
import dh_certified_tau as DH
from dh_certified_tau import build_tau_general, margin, dh_prime_side
c,N,prec=100,64,800
dw,dpos=dh_prime_side(c,prec)
print(f"Fine n_neg ramp around derived (k+J)_DH=0.61984 (c={c},N={N}); my const error ~5e-5:")
for kj in [0.610,0.614,0.616,0.618,0.6193,0.6198,0.6203,0.622,0.626]:
    A,DIM=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=dw,prime_positions=dpos)
    nn,pd,mg=margin(A,DIM)
    tag=" <-- derived +/-5e-4 band" if abs(kj-0.61984)<6e-4 else ""
    print(f"  (k+J)={kj:.4f}: n_neg={nn}{tag}",flush=True)

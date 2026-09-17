import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import ctx
from dh_exact import kj_exact
from dh_certified_tau import build_tau_general, dh_prime_side
from arb_ldlt_margin import inertia_with_margin
prec=800
for c in [25,50,100]:
    ctx.prec=prec
    kjd=float(kj_exact('0.75',5,c,prec).mid())
    dw,dpos=dh_prime_side(c,prec)
    A,DIM=build_tau_general(c,64,prec,nu=0.75,include_pole=False,kj_const=kjd,prime_weights=dw,prime_positions=dpos)
    npos,nneg,undet,mlow,mmid=inertia_with_margin(A,DIM)
    print(f"c={c}: n_pos={npos} n_neg={nneg} undetermined_pivot={undet} (None=fully certified) min|pivot|={float(mlow) if mlow else None}",flush=True)

import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_exact import kj_exact
from dh_certified_tau import build_tau_general
from arb_ldlt_margin import inertia_with_margin
def chi4(n):
    r=n%4; return 1 if r==1 else (-1 if r==3 else 0)
def side(c,prec):
    ctx.prec=prec; w=[];pos=[]
    for (q,p) in prime_powers_up_to(c):
        cq=chi4(q)
        if cq==0: continue
        w.append(arb(cq)*arb(p).log()*(arb(q)**arb('-0.5'))); pos.append(arb(q).log())
    return w,pos
for c,prec in [(50,800),(50,2500),(50,4000)]:
    ctx.prec=prec
    kj=float(kj_exact('0.75',4,c,prec).mid())
    w,pos=side(c,prec)
    A,DIM=build_tau_general(c,64,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=w,prime_positions=pos)
    npos,nneg,undet,mlow,mmid=inertia_with_margin(A,DIM)
    print(f"c={c} prec={prec}: n_neg={nneg} undetermined_pivot={undet} (None=certified) min|piv|={float(mlow) if mlow else None}",flush=True)

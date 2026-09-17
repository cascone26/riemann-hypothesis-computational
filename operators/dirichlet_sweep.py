import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_exact import kj_exact
from dh_certified_tau import build_tau_general, margin
def chi4(n):
    r=n%4; return 1 if r==1 else (-1 if r==3 else 0)
def side(c,prec):
    ctx.prec=prec; w=[];pos=[]
    for (q,p) in prime_powers_up_to(c):
        cq=chi4(q)
        if cq==0: continue
        w.append(arb(cq)*arb(p).log()*(arb(q)**arb('-0.5'))); pos.append(arb(q).log())
    return w,pos
c,prec=50,1000; ctx.prec=prec
w,pos=side(c,prec)
import math
print(f"L(chi4) c={c}: sweep arch constant (k+J) -> where is n_neg=0 (GRH-consistent)?")
print(f"  my kj_exact(nu=3/4,cond=4) = {float(kj_exact('0.75',4,c,prec).mid()):.5f}")
print(f"  kj_exact(nu=3/4,cond=16)   = {float(kj_exact('0.75',16,c,prec).mid()):.5f}  (cond=q^2?)")
print(f"  kj_exact(nu=3/4,cond=2)    = {float(kj_exact('0.75',2,c,prec).mid()):.5f}")
for kj in [0.70,0.75,0.80,0.8336,0.90,1.0,1.15,1.30,1.45]:
    A,DIM=build_tau_general(c,64,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=w,prime_positions=pos)
    nn,_,_=margin(A,DIM)
    print(f"    (k+J)={kj:.4f}: n_neg={nn}",flush=True)

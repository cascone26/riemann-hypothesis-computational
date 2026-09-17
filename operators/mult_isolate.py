import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general, margin
def chi4(n):
    r=n%4; return 1 if r==1 else (-1 if r==3 else 0)
def side(c,prec,flip_prime=None):
    ctx.prec=prec; w=[];pos=[]
    for (q,p) in prime_powers_up_to(c):
        cq=chi4(q)
        if cq==0: continue
        s=cq
        if flip_prime is not None and p==flip_prime:
            s=-cq   # break multiplicativity: flip this prime's whole tower sign
        w.append(arb(s)*arb(p).log()*(arb(q)**arb('-0.5'))); pos.append(arb(q).log())
    return w,pos
c,prec,N=50,1000,64; ctx.prec=prec
kj=0.90   # SAFE in the n_neg=0 plateau for the multiplicative case? test both at same arch
# use a constant where multiplicative L(chi4) is PSD: from sweep, need cond>=4.05; kj_exact(3/4,4.1)
from dh_exact import kj_exact
kj=float(kj_exact('0.75',4.1,c,prec).mid())
print(f"Same-arch (kj={kj:.4f}, cond=4.1) multiplicativity isolation at c={c}:")
for lbl,fp in [("multiplicative chi4",None),("break mult @ p=3",3),("break mult @ p=5",5),("break mult @ p=7",7)]:
    w,pos=side(c,prec,fp)
    A,DIM=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=w,prime_positions=pos)
    nn,_,_=margin(A,DIM)
    print(f"  {lbl:22s}: n_neg={nn}",flush=True)
print("If multiplicative=0 but breaking ONE prime's tower-sign -> n_neg>0, multiplicativity is load-bearing (fixed arch).")

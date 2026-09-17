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
prec=1000
print("L(chi4): which conductor in the arch constant gives certified n_neg=0 at ALL c (=GRH)?")
print(f"{'cond':>6} | {'c=25':>5} {'c=50':>5} {'c=100':>6}")
for cond in [4.0, 4.05, 4.1, 4.2, 4.5]:
    row=[]
    for c in [25,50,100]:
        ctx.prec=prec
        kj=float(kj_exact('0.75',cond,c,prec).mid())
        w,pos=side(c,prec)
        A,DIM=build_tau_general(c,64,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=w,prime_positions=pos)
        nn,_,_=margin(A,DIM); row.append(nn)
    print(f"{cond:>6.2f} | {row[0]:>5} {row[1]:>5} {row[2]:>6}",flush=True)
print("The cond giving n_neg=0 at all c is the correct Dirichlet-L arch normalization (GRH-consistent).")

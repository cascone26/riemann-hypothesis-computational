import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general, build_arb_tau, margin

def kmax_of(q,p):  # exponent a with q=p^a
    a=0;m=q
    while m>1: m//=p;a+=1
    return a

def zeta_side(c,prec,kmin=1,kmax=99):
    ctx.prec=prec
    pdata=[(q,p) for (q,p) in prime_powers_up_to(c) if kmin<=kmax_of(q,p)<=kmax]
    w=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    pos=[arb(q).log() for (q,p) in pdata]
    return w,pos

prec=800;N=64
print(f"CERTIFIED n_neg of zeta Weil form under Euler-product truncation in k (N={N}):")
print(f"{'c':>4} | {'FULL(all k)':>11} | {'k=1 only':>9} | {'k<=2':>6} | {'k<=3':>6} | {'k=2 only':>9}")
for c in [25,50,100]:
    ctx.prec=prec
    nf,_,_=margin(*build_arb_tau(c,N,prec))
    def run(kmin,kmax):
        w,pos=zeta_side(c,prec,kmin,kmax)
        A,D=build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=w,prime_positions=pos)
        return margin(A,D)[0]
    n1=run(1,1); n2=run(1,2); n3=run(1,3); nk2=run(2,2)
    print(f"{c:>4} | {nf:>11} | {n1:>9} | {n2:>6} | {n3:>6} | {nk2:>9}",flush=True)
print("If k<=2 ~ FULL(0) => the SQUARES (k=2) carry essentially all the prime-power slack.")

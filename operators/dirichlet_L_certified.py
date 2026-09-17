import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_exact import kj_exact
from dh_certified_tau import build_tau_general, build_arb_tau, margin
# L(s, chi_4): chi mod 4 (odd primitive), Euler product, conductor 4, NO pole, GRH holds.
# prime weights = chi(n) Lambda(n)/sqrt(n) on prime powers -- SIGNED (by chi), but multiplicative (Euler product).
def chi4(n):
    r=n%4
    return 1 if r==1 else (-1 if r==3 else 0)
def Lchi_prime_side(c,prec):
    ctx.prec=prec
    pdata=prime_powers_up_to(c)  # (q=p^k, p)
    w=[]; pos=[]
    for (q,p) in pdata:
        cq=chi4(q)                       # chi(p^k) = chi(p)^k (completely multiplicative)
        if cq==0: continue               # p=2 drops
        w.append(arb(cq)*arb(p).log()*(arb(q)**arb('-0.5')))
        pos.append(arb(q).log())
    return w,pos
prec,N=800,64
print(f"CERTIFIED n_neg for L(s,chi_4) (Euler product, SIGNED weights chi*Lambda, GRH holds) vs zeta, N={N}:")
print(f"{'c':>4} | {'zeta n_neg':>10} | {'L(chi4) n_neg':>13} | {'(k+J) used':>12}")
for c in [25,50,100]:
    ctx.prec=prec
    nz,_,_=margin(*build_arb_tau(c,N,prec))
    kj=float(kj_exact('0.75',4,c,prec).mid())   # nu=3/4 (odd chi), conductor 4
    w,pos=Lchi_prime_side(c,prec)
    A,DIM=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=w,prime_positions=pos)
    nn,pd,undet=margin(A,DIM)
    print(f"{c:>4} | {nz:>10} | {nn:>13} | {kj:>12.6f}",flush=True)
print("PREDICTION: L(chi4) n_neg=0 (GRH). If 0 -> EP forces PSD even with SIGNED weights => multiplicativity, not")
print("weight-sign, is the load-bearing ingredient (sharpens cycle 199).")

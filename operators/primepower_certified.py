import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
import dh_certified_tau as DH
from dh_certified_tau import build_tau_general, build_arb_tau, margin

def zeta_prime_side(c, prec, powers=True):
    ctx.prec=prec
    pdata=prime_powers_up_to(c)   # (q,p) with q=p^a<=c
    if not powers:
        pdata=[(q,p) for (q,p) in pdata if q==p]   # PRIMES ONLY (drop k>=2)
    w=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    pos=[arb(q).log() for (q,p) in pdata]
    return w,pos

prec=800
print("CERTIFIED test: does removing prime POWERS (k>=2) break zeta's Weil positivity?")
for c in [25,50,100]:
    ctx.prec=prec
    # full zeta (canonical) -- certified PSD baseline
    Af,Df=build_arb_tau(c,64,prec); nf,pdf,mf=margin(Af,Df)
    # zeta with primes ONLY
    w,pos=zeta_prime_side(c,prec,powers=False)
    Ap,Dp=build_tau_general(c,64,prec,nu=0.25,include_pole=True,prime_weights=w,prime_positions=pos)
    npr,pdp,mp_=margin(Ap,Dp)
    print(f"  c={c}: zeta FULL n_neg={nf} (cert PSD={pdf}) | zeta PRIMES-ONLY n_neg={npr} (PSD={pdp})",flush=True)
print("If PRIMES-ONLY n_neg>=1 while FULL=0 => prime powers (k>=2) are certifiably essential for positivity.")

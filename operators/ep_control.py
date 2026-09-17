import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general, build_arb_tau, margin
prec,N=800,64
for c in [25,100]:
    ctx.prec=prec
    nf,pf,_=margin(*build_arb_tau(c,N,prec))
    # reconstruct ALL prime powers via build_tau_general
    pdata=prime_powers_up_to(c)
    w=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    pos=[arb(q).log() for (q,p) in pdata]
    ng,pg,_=margin(*(build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=w,prime_positions=pos),)[0] if False else build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=w,prime_positions=pos))
    print(f"c={c}: build_arb_tau FULL n_neg={nf} | build_tau_general(all prime powers) n_neg={ng}  MATCH={nf==ng}",flush=True)

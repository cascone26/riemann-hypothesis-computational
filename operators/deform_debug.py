import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general, margin
def vonmangoldt(N):
    lam=[0.0]*(N+1)
    for p in range(2,N+1):
        if lam[p]==0:
            lp=math.log(p);pk=p
            while pk<=N: lam[pk]=lp;pk*=p
    return lam
c,N,prec=50,64,800; ctx.prec=prec
VM=vonmangoldt(c)
# Path A: prime-powers-only (known good, ep_control gave 0)
pdA=prime_powers_up_to(c)
wA=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdA]
posA=[arb(q).log() for (q,p) in pdA]
nA,_,_=margin(*build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=wA,prime_positions=posA))
# Path B: all-n with VM weights (deform t=0)
wB=[arb(VM[n])*(arb(n)**arb('-0.5')) for n in range(2,c+1)]
posB=[arb(n).log() for n in range(2,c+1)]
nB,_,_=margin(*build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=wB,prime_positions=posB))
print(f"Path A (prime-powers-only): n_neg={nA}")
print(f"Path B (all-n, VM weights):  n_neg={nB}")
# compare: are the nonzero weights the same set?
nzB=[(n,float(VM[n])) for n in range(2,c+1) if VM[n]>0]
print(f"Path B nonzero-weight n: {[n for n,_ in nzB][:20]}")
print(f"Path A q values:         {sorted([q for q,p in pdA])[:20]}")
print(f"VM check: VM[2]={VM[2]:.3f} VM[4]={VM[4]:.3f} VM[8]={VM[8]:.3f} VM[6]={VM[6]:.3f} VM[9]={VM[9]:.3f} VM[49]={VM[49]:.3f}")

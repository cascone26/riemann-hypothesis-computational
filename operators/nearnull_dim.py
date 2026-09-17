import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general, margin
from _dh_negative_control import lambda_D
N,prec=64,800
print(f"Near-null DIMENSION of zeta's Weil form (deformation jump at t=1e-6) vs prime cutoff c, N={N}:")
print(f"{'c':>4} | {'zeta n_neg(t=0)':>15} | {'near-null dim (n_neg @t=1e-6)':>28}")
for c in [25,50,100]:
    ctx.prec=prec
    ppset={q:p for (q,p) in prime_powers_up_to(c)}
    LD=lambda_D(c); pos=[arb(n).log() for n in range(2,c+1)]
    def wt(t): return [(arb(1-t)*(arb(ppset[n]).log() if n in ppset else arb(0)) + arb(t*float(LD[n])))*(arb(n)**arb('-0.5')) for n in range(2,c+1)]
    n0=margin(*build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=wt(0.0),prime_positions=pos))[0]
    ne=margin(*build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=wt(1e-6),prime_positions=pos))[0]
    print(f"{c:>4} | {n0:>15} | {ne:>28}",flush=True)
print("If near-null dim grows with c ~ #zeta-zeros resolved => the near-null modes ARE the resolved on-line zeros.")

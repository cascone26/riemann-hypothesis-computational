import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general, margin
from _dh_negative_control import lambda_D
c,N,prec=50,64,800; ctx.prec=prec
# Arb-EXACT zeta von Mangoldt: Lambda(n)=arb(p).log() on prime powers, else 0
ppset={q:p for (q,p) in prime_powers_up_to(c)}
def Lam_zeta_arb(n):
    return arb(ppset[n]).log() if n in ppset else arb(0)
LD=lambda_D(c)
pos=[arb(n).log() for n in range(2,c+1)]
print(f"Arb-EXACT zeta weights; deform w(t)=(1-t)Lambda_arb + t*Lambda_D, all n<=c={c}, N={N}")
print(f"{'t':>7} | {'n_neg':>6}")
for t in [0.0,1e-6,1e-4,0.001,0.01,0.05,0.1,0.25,0.5,1.0]:
    w=[(arb(1-t)*Lam_zeta_arb(n) + arb(t*float(LD[n])))*(arb(n)**arb('-0.5')) for n in range(2,c+1)]
    A,D=build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=w,prime_positions=pos)
    nn,_,_=margin(A,D)
    print(f"{t:>7.0e} | {nn:>6}" if t<0.01 else f"{t:>7.2f} | {nn:>6}",flush=True)
print("t=0 must be 0 (Arb-exact zeta control). Smallest t with n_neg>0 = zeta's PSD-margin width toward non-mult.")

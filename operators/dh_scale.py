import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import ctx
import dh_certified_tau as DH
from dh_certified_tau import build_tau_general, build_arb_tau, margin, dh_prime_side
c,prec=14,800
print(f"SCALING: does D-H n_neg grow robustly (const-insensitive) while zeta stays 0?  c={c}")
print(f"{'N':>4} | {'zeta n_neg(exact)':>17} | {'D-H n_neg @0.5945':>17} | {'D-H n_neg @0.5957':>17}")
for N in [32,48,64,96]:
    dw,dpos=dh_prime_side(c,prec)
    Az,Dz=build_arb_tau(c,N,prec); nz,_,_=margin(Az,Dz)
    Aa,Da=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=0.5945,prime_weights=dw,prime_positions=dpos); na,_,_=margin(Aa,Da)
    Ab,Db=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=0.5957,prime_weights=dw,prime_positions=dpos); nb,_,_=margin(Ab,Db)
    print(f"{N:>4} | {nz:>17} | {na:>17} | {nb:>17}",flush=True)
print("robust discrimination = D-H n_neg grows with N AND stays >=1 at BOTH constants, while zeta stays 0.")

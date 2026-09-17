# NEW bridge: does the pole-less Weil form Q0 have EXACTLY ONE negative eigenvalue for MULTIPLICATIVE
# L-functions (kappa=1, pole lifts it -> PSD) but MORE for NON-mult (D-H's off-line pairs survive)?
import sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
sys.path.insert(0,'/Users/scones/projects/riemann/external/groskin_2607.02828/code')
from flint import ctx
from arb_ldlt_certify import build_arb_tau
from dh_certified_tau import build_tau_general, dh_prime_side, margin
from dh_exact import kj_exact
print(f"{'c':>4} {'fn':>5} {'n_neg(Q0,no pole)':>17} {'n_neg(full,+pole)':>17} {'lifted by pole':>14}")
for c in [50,100]:
    N=64;prec=800;ctx.prec=prec
    # zeta: nu=1/4. no-pole vs with-pole
    Az0,Dz0=build_tau_general(c,N,prec,nu=0.25,include_pole=False)
    Azf,Dzf=build_tau_general(c,N,prec,nu=0.25,include_pole=True)
    nz0=margin(Az0,Dz0)[0]; nzf=margin(Azf,Dzf)[0]
    print(f"{c:>4} {'zeta':>5} {nz0:>17} {nzf:>17} {nz0-nzf:>14}",flush=True)
    # D-H: nu=3/4, cond=5, signed Lambda_D. no-pole (D-H is entire, but we ADD the same pole term to test lift)
    kj=float(kj_exact('0.75',5,c,prec).mid()); dw,dpos=dh_prime_side(c,prec)
    Ad0,Dd0=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=dw,prime_positions=dpos)
    Adf,Ddf=build_tau_general(c,N,prec,nu=0.75,include_pole=True,kj_const=kj,prime_weights=dw,prime_positions=dpos)
    nd0=margin(Ad0,Dd0)[0]; ndf=margin(Adf,Ddf)[0]
    print(f"{c:>4} {'D-H':>5} {nd0:>17} {ndf:>17} {nd0-ndf:>14}",flush=True)
print("Bridge: if zeta Q0 n_neg=1 -> full=0 (pole lifts the ONE), and D-H Q0 n_neg>1 -> full>0 (pole lifts")
print("only one, off-line pairs survive), then MULTIPLICATIVITY => kappa=1 => pole suffices. Ties discriminator to pole-lift.")

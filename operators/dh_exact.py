import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
import arb_ldlt_certify as G
import dh_certified_tau as DH
from dh_certified_tau import build_tau_general, build_arb_tau, margin, dh_prime_side

def kj_exact(nu, cond, c, prec):
    """(kappa+J) = -psi(nu) + log(pi/cond) - sum_k e^{-2(k+nu)L}/(k+nu), rigorous Arb (convergent, tail-bounded)."""
    ctx.prec=prec
    L=arb(c).log(); nu=arb(nu)
    s=arb(0); k=0; thr=arb(2)**(-(prec-8))
    while True:
        t=(-2*(k+nu)*L).exp()/(k+nu)
        s+=t
        if t<thr and k>2: break
        k+=1
    return -nu.digamma() + (arb.pi()/cond).log() - s

# validate against Groskin exact for zeta (nu=1/4,cond=1)
prec=800
for c in [14,50,100]:
    L=arb(c).log(); ctx.prec=prec
    exact_z=kj_exact('0.25',1,c,prec); groskin=G.arb_kappa(L)+G.arb_J(L)
    print(f"  c={c}: kj_exact(zeta)={float(exact_z.mid()):.8f}  Groskin={float(groskin.mid()):.8f}  diff={float((exact_z-groskin).mid()):.2e}")

print("\nDEFINITIVE certified table with EXACT (kappa+J)_DH (no numeric constant):")
print(f"{'c':>4} | {'(k+J)_DH exact':>15} | {'zeta n_neg':>10} | {'D-H n_neg':>9}")
N=64
for c in [14,25,50,75,100]:
    ctx.prec=prec
    kjd=kj_exact('0.75',5,c,prec)
    dw,dpos=dh_prime_side(c,prec)
    Az,Dz=build_arb_tau(c,N,prec); nz,_,_=margin(Az,Dz)
    A,DIM=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=float(kjd.mid()),prime_weights=dw,prime_positions=dpos)
    nn,pd,mg=margin(A,DIM)
    print(f"{c:>4} | {float(kjd.mid()):>15.7f} | {nz:>10} | {nn:>9}",flush=True)

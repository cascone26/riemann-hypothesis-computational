# Rigorous Arb certification of the genuine E11a (weight-2 level-11 elliptic-curve) Weil form (cycle 268).
# Exact integer a_p -> Satake b_{p^k} in Arb; arch via acb.integral; Omega=Re psi(1+ir)-log2pi+0.5 log11.
import sys, math
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
sys.path.insert(0,'/Users/scones/projects/riemann/external/groskin_2607.02828/code')
import numpy as np
from flint import arb, acb, arb_mat, ctx
from dh_certified_tau import margin
from sympy import primerange
ctx.prec=300
PI=arb.pi(); SQRT2PI=(2*PI).sqrt(); LOG2PI=(2*PI).log(); HALF_LOG11=arb(11).log()/2
atil=np.load('E11_an.npy')
def ap_int(p): return round(float(atil[p])*math.sqrt(p))
NM=2500
# prime side: Satake power sums b_{p^k}=alpha^k+beta^k, alpha+beta=atil_p=a_p/sqrt(p), alpha*beta=1 (good);
# bad prime 11: b_{11^k}=atil_11^k (single Satake). Lambda_E(p^k)=b_{p^k} log p.
PRIMES=[]
for p in primerange(2,NM):
    ap=arb(ap_int(p)); atp=ap/arb(p).sqrt(); lp=arb(p).log()
    if p==11:
        bcur=arb(1)  # atil_11^0=1; we build b_{11^k}=atil_11^k iteratively
        k=1; pk=p
        while pk<NM:
            bcur=bcur*atp  # atil_11^k
            PRIMES.append((pk,lp,bcur)); k+=1; pk*=p
    else:
        bprev,bcur=arb(2),atp; k=1; pk=p
        while pk<NM:
            if k>=2:
                bnew=atp*bcur-bprev; bprev,bcur=bcur,bnew
            PRIMES.append((pk,lp,bcur)); k+=1; pk*=p
def g_arb(aj,ak,y):
    s=aj+ak; d=aj-ak; sp2=PI.sqrt()/2
    return (-y*y/4).exp()*sp2*((-s*s/4).exp()*(d*y/2).cos()+(-d*d/4).exp()*(s*y/2).cos())
def prime_arb(aj,ak):
    tot=arb(0)
    for (n,lp,b) in PRIMES:
        y=arb(n).log()
        tot+=b*lp*(arb(n)**arb('-0.5'))*g_arb(arb(aj),arb(ak),y)
    tail=arb(0).union(arb('3e-5'))  # NM=2500 (E11_an cap); tail VERIFIED <1.12e-5 (x2:2.23e-5), +-1.5e-5 covers
    return tot+tail-arb('3e-5')/2
def phihat(a,r):
    A=acb(a); return (SQRT2PI/2)*((-(r-A)*(r-A)/2).exp()+(-(r+A)*(r+A)/2).exp())
def arch_arb(aj,ak):
    R=arb(60)
    def integrand(r,_):
        z=acb(1)+acb(0,1)*r
        Om=z.digamma().real - acb(LOG2PI) + acb(HALF_LOG11)
        return phihat(aj,r)*phihat(ak,r)*Om
    return (acb.integral(integrand, acb(-R), acb(R)).real)/(2*PI)
def Q(aj,ak): return 2*(arch_arb(aj,ak)-prime_arb(aj,ak))
A=[8.0,12.0,16.0,20.0,24.0]
n=len(A); M=arb_mat(n,n)
print(f"certified E11a Weil form, prec={ctx.prec}, NM={NM}...",flush=True)
for i in range(n):
    for j in range(i,n):
        q=Q(A[i],A[j]); M[i,j]=q; M[j,i]=q
        print(f"  Q[{A[i]:.0f},{A[j]:.0f}]={float(q.mid()):+.5f} (rad {float(q.rad()):.1e})",flush=True)
nn,pd,mp=margin(M,n)
print(f"\nCERTIFIED E11a Weil form: n_neg={nn} PD={pd} min_pivot={float(mp) if mp is not None else None}",flush=True)

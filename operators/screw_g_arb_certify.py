import sys; sys.path.insert(0,"/Users/scones/projects/riemann/operators")
from flint import arb, ctx
import numpy as np, mpmath as mp
ctx.prec=200; mp.mp.dps=55
def A(x): return arb(str(x))                       # mpmath -> arb (tight)
LOGPI=A(mp.log(mp.pi)); PSI_Q=A(mp.digamma(mp.mpf(1)/4))
C_ARCH=PSI_Q-LOGPI
PHI1=A(mp.polygamma(1,mp.mpf(1)/4))                 # trigamma(1/4)=Phi(1,2,1/4)
def Phi_2_quarter(z):                               # sum z^n/(n+1/4)^2, rigorous geometric tail
    tot=arb(0); zn=arb(1); n=0; zf=float(z)
    while True:
        tot+=zn/(arb(n)+A('0.25'))**2
        if n>3 and zf<1:
            tb=(zf**(n+1))/((n+1.25)**2)/(1-zf)
            if tb<1e-58: tot+=arb(0).union(arb(tb)); break
        zn*=z; n+=1
        if n>500000: break
    return tot
def vonmangoldt(N):
    Lam=np.zeros(N+1); comp=np.zeros(N+1,bool)
    for p in range(2,N+1):
        if not comp[p]:
            comp[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: Lam[pk]=lp; pk*=p
    return Lam
def g_arb(t):
    T=A(t)
    boundary=-4*((T/2).exp()+(-T/2).exp()-2)
    arch_core=-(T/2)*C_ARCH
    arch_tail=-A('0.25')*(PHI1-(-T/2).exp()*Phi_2_quarter((-2*T).exp()))
    N=int(np.floor(np.exp(t))); Lam=vonmangoldt(N); prime=arb(0)
    for n in range(2,N+1):
        if Lam[n]>0:
            prime+=A(Lam[n])/arb(n).sqrt()*(T-arb(n).log())
    return boundary+prime+arch_core+arch_tail
print("RIGOROUS Arb g(t)<0 (partial Cor 4.3), prec=200:")
for t in [1,2,3,5,8,12]:
    g=g_arb(t); mid=float(g.mid()); rad=float(g.rad())
    print(f"  t={t:>2}: g = {mid:+.7f} +/- {rad:.1e}   {'CERTIFIED g<0' if (mid+rad)<0 else 'NOT cert'}",flush=True)

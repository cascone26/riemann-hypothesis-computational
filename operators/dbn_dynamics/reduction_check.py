#!/usr/bin/env python3
"""§199.3a(v2): validate the 1D-integral REDUCTION Q[k,l]=∫_0^{2a} g(t)D_kl(t)dt (D_kl=C_kl+C_lk, cos-overlap),
via fast mpmath quad (split at prime kink log2), vs mpmath weil_matrices. No sympy. If it matches, the reduction
that the rigorous Arb build rests on is correct."""
import sys, mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, screw_g_mp, prime_data
mp.mp.dps=30
a=mp.mpf('0.549'); pl=prime_data(a); log2=mp.log(2); b=mp.pi/(2*a)
def dprime(k,x): return (k*mp.pi/(2*a))*mp.cos(k*mp.pi*(x+a)/(2*a))
def C(k,l,t):   # ∫_{t-a}^a d_k(x) d_l(x-t) dx
    return mp.quad(lambda x: dprime(k,x)*dprime(l,x-t), [t-a, a])
def D(k,l,t): return C(k,l,t)+C(l,k,t)
def Q(k,l):
    f=lambda t: screw_g_mp(t,pl,a)[0]*D(k,l,t)
    return mp.quad(f,[0, float(log2), float(2*a)])
modes=[1,3,5,7]
Qf,Qa,Mass,Vv,xs,w=weil_matrices(a,24,1000)
ev=[c for c in range(24) if (c+1)%2==1][:4]
print("§199.3a(v2) 1D-reduction check vs mpmath @ a=0.549")
maxe=mp.mpf(0)
for ik,k in enumerate(modes):
    for il,l in enumerate(modes):
        if il<ik: continue
        val=Q(k,l); ref=Qf[ev[ik],ev[il]]; e=abs(val-ref); maxe=max(maxe,e)
        print(f"  Q[{k},{l}]  1Dquad={mp.nstr(val,12)}  weil_mtx={mp.nstr(ref,12)}  err={mp.nstr(e,3)}",flush=True)
print(f"\nmax err = {mp.nstr(maxe,3)} (tiny ⟹ reduction CORRECT; rigorous Arb version = same integral in acb.integral)")

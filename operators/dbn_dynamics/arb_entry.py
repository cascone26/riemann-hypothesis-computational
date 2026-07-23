#!/usr/bin/env python3
"""§199.2: rigorous Arb entry Q[k,l]=∫∫ g(x-y)d_k(x)d_l(y). Reduce to ∫_0^{2a} g(t)[C_kl(t)+C_lk(t)] dt (g even),
C_kl(t)=∫ d_k(x)d_l(x-t)dx (elementary, sympy). Integrate rigorously (acb.integral) split at the prime kink t=log2.
Validate (1,1) entry vs mpmath weil_matrices Qf[0,0]."""
import sys
import sympy as sp
from flint import acb, arb, ctx
ctx.prec=260
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
import mpmath as mp
from mpmath_weil_hiprec import weil_matrices, prime_data
from arb_screw import g_arb

# --- sympy: C_kl(t) for t in [0,2a], overlap x in [t-a, a] ---
xs,ts,as_,ks,ls = sp.symbols('x t a k l', real=True, positive=True)
dk=(ks*sp.pi/(2*as_))*sp.cos(ks*sp.pi*(xs+as_)/(2*as_))
dl=(ls*sp.pi/(2*as_))*sp.cos(ls*sp.pi*(xs-ts+as_)/(2*as_))
Ckl=sp.integrate(dk*dl,(xs,ts-as_,as_))          # closed form in t (params a,k,l)
Ckl=sp.simplify(Ckl)
def Cfun(k,l,aval):
    e=Ckl.subs({ks:k,ls:l,as_:sp.Rational(str(aval))})
    f=sp.lambdify(ts,e,modules=[{'sin':lambda z:z.sin(),'cos':lambda z:z.cos(),'pi':acb.pi()},'math'])
    return f

if __name__=="__main__":
    aval=0.549; mp.mp.dps=45; am=mp.mpf(str(aval))
    pl=prime_data(am); log2=float(mp.log(2)); twoa=2*aval
    k=l=1
    Ck=Cfun(k,l,aval); Cl=Cfun(l,k,aval)
    def integrand(t,_):
        gt=g_arb(float(t.real.mid()) if hasattr(t,'real') else float(t), pl, am)  # g at acb point (real)
        return gt*(Ck(t)+Cl(t))
    # split at kink t=log2
    I1=acb.integral(integrand, acb(0), acb(log2))
    I2=acb.integral(integrand, acb(log2), acb(twoa))
    Q11_arb=I1+I2
    Qf,Qa,Mass,Vv,xs2,w=weil_matrices(am,24,1200)
    print(f"§199.2 rigorous (1,1) entry @ a={aval}")
    print(f"  mpmath Qf[0,0] = {mp.nstr(Qf[0,0],18)}")
    print(f"  Arb rigorous  = {str(Q11_arb.real)}")
    lo=float(Q11_arb.real.mid())-float(Q11_arb.real.rad()); hi=float(Q11_arb.real.mid())+float(Q11_arb.real.rad())
    print(f"  contains mpmath? {lo <= float(Qf[0,0]) <= hi}")

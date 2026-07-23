#!/usr/bin/env python3
"""§199 increment 1: the g_a screw kernel g(t) in rigorous Arb intervals. Validate against mpmath_weil_hiprec.screw_g_mp.
g(t)=boundary -4(e^{t/2}+e^{-t/2}-2) + arch_core -(|t|/2)(ψ(1/4)-logπ) + arch_tail -(1/4)[Φ0 - e^{-|t|/2}Φ(e^{-2|t|})]
      + prime Σ_{logn≤|t|} Λ(n)/√n (|t|-logn).   Φ(z)=Σ z^n/(n+1/4)^2 (Lerch), Φ0=ζ(2,1/4)."""
import sys
from flint import arb, acb, ctx
ctx.prec = 300
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import screw_g_mp, prime_data

def hurwitz2(x):  # ζ(2,x)=Σ_{n≥0}(n+x)^-2, rigorous via acb Hurwitz zeta
    return acb(2).zeta(acb(x)).real

def lerch2(z, x, terms):  # Σ_{n≥0} z^n/(n+x)^2, z=e^{-2|t|}<1, rigorous tail bound
    s = arb(0); zp = arb(1)
    for n in range(terms):
        s += zp/(arb(n)+x)**2; zp *= z
    # tail: Σ_{n≥terms} z^n/(n+x)^2 ≤ z^terms/(terms+x)^2 * 1/(1-z)  (enclose)
    tail = zp/(arb(terms)+x)**2/(arb(1)-z)
    return s + arb(0).union(tail)   # enclose [s, s+tail]

def g_arb(t, primes, a):
    at = abs(arb(t)); half=arb(1)/2; quarter=arb(1)/4
    boundary = -4*((at*half).exp()+(-at*half).exp()-2)   # even in t; use |t| ok since even
    om0 = quarter.digamma() - arb.pi().log()
    arch_core = -(at/2)*om0
    phi0 = hurwitz2(quarter)
    z = (-2*at).exp()
    phit = lerch2(z, quarter, 200)
    arch_tail = -(quarter)*(phi0 - (-at/2).exp()*phit)
    prime = arb(0)
    for (Lam, logn, sqrtn) in primes:
        if float(logn) <= float(at.mid()):   # logn are exact-ish; guard
            prime += arb(str(Lam))/arb(str(sqrtn))*(at-arb(str(logn)))
    return boundary+arch_core+arch_tail+prime

if __name__=="__main__":
    mp.mp.dps=40; a=mp.mpf('0.549')
    pl=prime_data(a)
    parb=[(str(L),str(ln),str(sn)) for (L,ln,sn) in pl]
    print("§199.1 Arb screw kernel validation vs mpmath (should agree; Arb gives an interval):")
    print(f"{'t':>7} {'mpmath g(t)':>22} {'Arb g(t) (interval)':>32}")
    for tt in [0.1,0.5,0.7,0.9]:
        gm,_=screw_g_mp(mp.mpf(str(tt)),pl,a)
        ga=g_arb(tt,pl,a)
        print(f"{tt:>7} {mp.nstr(gm,16):>22} {str(ga):>32}")
    print("\nREAD: Arb interval must CONTAIN the mpmath value at each t ⟹ kernel port correct. Then build rigorous entries.")

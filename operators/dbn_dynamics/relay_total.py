#!/usr/bin/env python3
"""§190: SHARP relay measure — TOTAL arch negativity Σ|neg eigenvalues|(a) vs budget B(a)=ΣΛ(n)/√n.
§189 showed |α₁|/B is a sawtooth peaking ~0.85 (single deepest mode = necessary-not-sufficient). Test the
total negativity (all arch negative modes summed): does Σ|neg λ|/B ride ~1 with the next prime entering as it→1?
Also #neg(arch)(a). If tight ⟹ the relay Σ(prime norms) ≈ Σ(arch negativity) is the exact anatomy of A+B≥P."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
from relay_test import prime_budget
mp.mp.dps=40; K=24; M=600
avals=[mp.mpf(x) for x in '0.45,0.549,0.60,0.693,0.75,0.805,0.85'.split(',')]
print(f"§190 total-negativity relay | dps=40 K={K} M={M}")
print(f"{'a':>6} {'#neg':>5} {'Σ|neg λ|(arch)':>15} {'budget B':>11} {'Σ|neg|/B':>10} {'full λmin':>12}")
for a in avals:
    Qf,Qa,Mass,Vv,xs,ww=weil_matrices(a,K,M)
    Ea,_=gen_eig_sym(Qa,Mass); Ef,_=gen_eig_sym(Qf,Mass)
    negs=[e for e in Ea if e<0]; tot=mp.fsum(-e for e in negs)
    B,_=prime_budget(a); r=tot/B if B>0 else mp.mpf('nan')
    print(f"{mp.nstr(a,4):>6} {len(negs):>5} {mp.nstr(tot,7):>15} {mp.nstr(B,6):>11} {mp.nstr(r,5):>10} {mp.nstr(Ef[0],4):>12}",flush=True)
print("\nREAD: Σ|neg λ|/B ~1 & →1 at each prime entry ⟹ EXACT relay (Σ prime norms = Σ arch negativity). Wanders ⟹ alignment-limited.")

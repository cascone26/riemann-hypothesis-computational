#!/usr/bin/env python3
"""§201c(A): confirm D_kl(t) -> -2kβ·sin(lβt) as l->inf (k fixed). Pure closed-form (Dcf), fast."""
import sys; sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
import mpmath as mp
from dkl_closedform import Dcf, a, b
mp.mp.dps=30
print(f"a={mp.nstr(a,4)} b=β={mp.nstr(b,5)}   D_kl(t) / (-2kβ sin(lβt)) -> 1 ?")
for k in [1,3]:
    for t in [mp.mpf('0.3'),mp.mpf('0.7')]:
        print(f" k={k} t={mp.nstr(t,2)}: " + "  ".join(
            f"l={l}:{mp.nstr(Dcf(k,l,t)/(-2*k*b*mp.sin(l*b*t)),6)}" for l in [21,41,81,161,321]))
print("READ: ratios -> 1 ⟹ D_kl ~ -2kβ sin(lβt) exact (backs the 1/l coupling bound C_k=2k[g(0)+g(2a)]).")

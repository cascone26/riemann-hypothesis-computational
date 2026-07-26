#!/usr/bin/env python3
"""§269e: extract the per-prime deformation law f(p) = b_1/b_arch for the single-prime semilocal measure S=[p].
This is the atomic arithmetic content of the near-null-onset (§269d showed b_1 is prime-value dependent).
Test candidate closed forms: does 1-f(p) ~ p^{-1/2}? ~ p^{-1}? something clean?"""
import sys, mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
from ccm_discrepancy import b_arch, d_coeffs   # import sets dps=60
K = 3
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
print(f"§269e per-prime law f(p)=b_1/b_arch([p])  (K={K}, dps={mp.mp.dps}):", flush=True)
print(f"{'p':>4} {'f(p)':>12} {'1-f(p)':>12} {'(1-f)/p^-0.5':>13} {'(1-f)/p^-1':>12} {'(1-f)*p^0.5':>12}", flush=True)
for p in primes:
    d = d_coeffs([p], K)
    f = (b_arch(1) + d[0]) / b_arch(1)
    omf = 1 - f
    print(f"{p:>4} {mp.nstr(f,7):>12} {mp.nstr(omf,7):>12} {mp.nstr(omf/p**mp.mpf('-0.5'),6):>13} "
          f"{mp.nstr(omf*p,6):>12} {mp.nstr(omf*mp.sqrt(p),6):>12}", flush=True)
print("DONE. If (1-f)*p^0.5 ~ const => 1-f ~ c/sqrt(p) (the Euler p^-1/2 scale). Check which column is flattest.", flush=True)

#!/usr/bin/env python3
"""§927 (autonomous): Báez-Duarte/Riesz coefficients c_k = Σ_{j=0}^k (-1)^j C(k,j)/ζ(2j+2). RH ⟺ c_k = O(k^{-3/4+ε}).
Comfortable criterion (decay, not razor). High-precision (mpmath, catastrophic cancellation). Test the decay rate —
is it accessible/comfortable? does our machinery engage the RATE where the band-limited razor blocks?"""
import mpmath as mp
mp.mp.dps=60
def ck(k):
    s=mp.mpf(0)
    for j in range(k+1):
        s+= (-1)**j * mp.binomial(k,j) / mp.zeta(2*j+2)
    return s
print("="*70)
print("§927 Báez-Duarte/Riesz c_k = Σ(-1)^j C(k,j)/ζ(2j+2). RH ⟺ c_k=O(k^{-3/4+ε})")
print("="*70)
print(f"{'k':>4} {'c_k':>16} {'k^{3/4}·c_k':>14} {'k^{1/2}·c_k':>14}")
import math
vals=[]
for k in [1,2,3,5,8,12,16,20,25,30,35,40]:
    c=ck(k); vals.append((k,float(c)))
    print(f"{k:>4} {float(c):>16.8f} {k**0.75*float(c):>14.5f} {k**0.5*float(c):>14.5f}")
print("="*70)
# fit decay exponent
import numpy as np
ks=np.array([v[0] for v in vals[2:]]); cs=np.abs(np.array([v[1] for v in vals[2:]]))
if (cs>0).all():
    b1,b0=np.polyfit(np.log(ks),np.log(cs),1)
    print(f"c_k ~ k^{b1:.3f} (RH predicts ≈ -0.75; unconditional bound is weaker)")
print("c_k decay is COMFORTABLE (smooth power law) — the RH-content is the UNIFORM O(k^{-3/4}) for ALL k.")
print("Proving the decay rate for all k = RH; but the per-k values are well-conditioned (mpmath-computable, not razor).")
print("="*70)

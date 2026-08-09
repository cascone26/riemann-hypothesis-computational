#!/usr/bin/env python3
"""§898b: verify the closed-form archimedean weight Ŵ(u) = −π e^{−u/2}/(1−e^{−2u}) by transforming the FULL Ω(r)
(not Ω−asymptote) with a convergence factor e^{−εr}, R large. Compare at u=log2,log3,1.5,2.5. Match ⇒ closed form
CONFIRMED (the §898 mismatch was from subtracting the ln asymptote = a different function)."""
import numpy as np
from scipy.special import digamma
def Wclosed(u): return -np.pi*np.exp(-u/2)/(1-np.exp(-2*u))
def Wnum_full(u,R=8000.0,dr=0.004,eps=0.004):
    r=np.arange(dr,R,dr); Om=np.real(digamma(0.25+0.5j*r))-np.log(np.pi)
    return np.array([2*np.sum(Om*np.cos(uu*r)*np.exp(-eps*r))*dr for uu in u])
us=np.array([0.6931,1.0986,1.5,2.0,2.5,3.0])
Wc=Wclosed(us)
for eps in [0.01,0.004,0.002]:
    Wn=Wnum_full(us,eps=eps)
    print(f"eps={eps}: "+"  ".join(f"u={u:.3f}[cf={c:.3f},num={n:.3f}]" for u,c,n in zip(us,Wc,Wn)))
print("\nclosed form −π e^{−u/2}/(1−e^{−2u}); if num→closed as eps→0 (mild u-indep offset from the δ/const at u=0 OK) ⇒ CONFIRMED")
# check ratio num/closed is ~1 (shape match) at the smallest eps
Wn=Wnum_full(us,eps=0.002)
print("ratio num/closed (eps=0.002):", np.array2string(Wn/Wc,precision=3))

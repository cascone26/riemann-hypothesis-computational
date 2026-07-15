#!/usr/bin/env python3
"""
beurling_selberg_test.py — Does a pointwise Selberg majorant close the wall A+B >= P?

Archimedean kernel on the log-line (Weil explicit formula, derived + verified):
  Omega_hat(x) = -e^{-|x|/2}/(1-e^{-2|x|})     [regular part; verified vs numeric inv-FT to 0.1%]
So the archimedean density magnitude at x=log k is  K(log k) = k^{3/2}/(k^2-1).
Prime comb weight at x=log k is  Lambda(k)/sqrt(k).

RESULT: ratio = [Lambda(k)/sqrt(k)] / [k^{3/2}/(k^2-1)] = Lambda(k)(1-1/k^2) -> log(p) -> inf.
Prime comb EXCEEDS archimedean kernel at every prime >= 5. No pointwise majorant exists.
Positivity of A+B-P is therefore GLOBAL (zeros constrain prime alignment), not pointwise.
"""
import numpy as np
from scipy.special import digamma
import mpmath as mp
LOGPI=np.log(np.pi)
def Omega(r): return np.real(digamma(0.25+0.5j*r))-LOGPI
def K_analytic(x): return np.exp(-x/2)/(1-np.exp(-2*x))
def trap(y,x): return float(np.sum((y[1:]+y[:-1])*0.5*np.diff(x)))
def Omega_hat_num(x,R=6000,npts=6000000):
    r=np.linspace(1e-6,R,npts); tail=np.log(r/(2*np.pi))
    return trap((Omega(r)-tail)*np.cos(r*x),r)/np.pi - 1/(2*x)
def vonmangoldt(k):
    n,p,f=k,0,2
    while f*f<=n:
        if n%f==0:
            p=f
            while n%f==0:n//=f
            if n>1: return 0.0
            break
        f+=1
    if n>1: p=n
    return float(mp.log(p)) if p>0 else 0.0

if __name__=="__main__":
    print("Kernel check (analytic vs numeric inv-FT):")
    for x in [0.5,0.693147,1.098612,1.609438]:
        print(f"  x={x:.3f}: num={Omega_hat_num(x):.5f}  -K={-K_analytic(x):.5f}")
    print("\nPrime comb vs archimedean kernel:")
    for k in range(2,30):
        L=vonmangoldt(k)
        if L==0: continue
        ratio=L*(1-1/k**2)
        print(f"  k={k:>3} ratio={ratio:.4f} {'PRIME EXCEEDS' if ratio>1 else 'arch>=prime'}")

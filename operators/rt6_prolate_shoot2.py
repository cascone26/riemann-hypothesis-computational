#!/usr/bin/env python3
"""§463 (faster shooting): φ(1;μ)=0 = regular-at-1 condition. Coarse scan + refine."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy'); FP=(2*np.pi)**2
def flux_at_1(E):
    mu=-E*E
    def rhs(x,y): return [y[1]/(1-x*x), (FP*x*x-mu)*y[0]]
    s=solve_ivp(rhs,[0,1-1e-5],[1.0,0.0],rtol=1e-9,atol=1e-11,max_step=0.02)
    return s.y[1,-1]
Es=np.linspace(3,35,320)
fp=np.array([flux_at_1(E) for E in Es])
roots=[]
for i in range(len(Es)-1):
    if fp[i]*fp[i+1]<0:
        a,b=Es[i],Es[i+1]
        for _ in range(40):
            m=(a+b)/2
            if flux_at_1(a)*flux_at_1(m)<0: b=m
            else: a=m
        roots.append((a+b)/2)
roots=np.array(roots)
print(f"§463 regular-at-1 E-values (φ(1)=0) vs zeta zeros:")
print(f"{'k':>3} {'E=√|μ|':>10} {'γ_k':>10} {'ratio':>7} {'E-spacing':>10}")
for k in range(min(18,len(roots))):
    gz=Z[k] if k<len(Z) else np.nan
    sp=roots[k]-roots[k-1] if k>0 else 0
    print(f"{k+1:>3} {roots[k]:10.4f} {gz:10.4f} {roots[k]/gz if gz==gz else float('nan'):7.4f} {sp:10.4f}")
print(f"\n#roots in [3,35]={len(roots)}, #zeros<35={int(np.sum(Z<35))}, #zeros in [3,35]... first zero=14.13")
print("reading: match to γ_k => correct realization; right count but shifted => need the ∞ B=0 condition too.")

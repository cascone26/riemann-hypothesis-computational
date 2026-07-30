#!/usr/bin/env python3
"""§463c: cleaner B-extraction. Integrate regular solution from x=1⁺ to X=70; extract B by LEAST-SQUARES fit of
u=x·ξ to A sin(2πx)+B cos(2πx) over the last 6 periods [64,70] (kills 1/x² contamination). Find μ<0 where B=0."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy'); FP=(2*np.pi)**2
def B_clean(E):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    e=1e-4
    xs=np.linspace(64,70,600)
    s=solve_ivp(rhs,[1+e,70],[1.0,(mu-FP)/2.0],rtol=1e-9,atol=1e-11,max_step=0.02,t_eval=xs)
    u=xs*s.y[0]
    Mm=np.column_stack([np.sin(2*np.pi*xs),np.cos(2*np.pi*xs)])
    (A,B),*_=np.linalg.lstsq(Mm,u,rcond=None)
    return B
Es=np.linspace(3,30,180)
Bs=np.array([B_clean(E) for E in Es])
roots=[]
for i in range(len(Es)-1):
    if Bs[i]*Bs[i+1]<0:
        a,b=Es[i],Es[i+1]
        for _ in range(30):
            m=(a+b)/2
            if B_clean(a)*B_clean(m)<0: b=m
            else: a=m
        roots.append((a+b)/2)
roots=np.array(roots)
print(f"§463c clean B-fit (X=70, 6-period lstsq). E=√|μ| vs zeros:")
print(f"{'k':>3} {'E':>9} {'γ_k':>9} {'ratio':>7}")
for k in range(min(15,len(roots))):
    gz=Z[k] if k<len(Z) else np.nan
    print(f"{k+1:>3} {roots[k]:9.4f} {gz:9.4f} {roots[k]/gz if gz==gz else float('nan'):7.4f}")
print(f"#roots [3,30]={len(roots)}, #zeros<30={int(np.sum(Z<30))}")

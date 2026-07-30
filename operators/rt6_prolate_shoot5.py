#!/usr/bin/env python3
"""§463c-fast: clean B-fit (X=45, 5-period lstsq), coarse scan, to settle extraction-vs-extension."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy'); FP=(2*np.pi)**2
def B_clean(E):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    xs=np.linspace(40,45,400)
    s=solve_ivp(rhs,[1+1e-4,45],[1.0,(mu-FP)/2.0],rtol=1e-8,atol=1e-10,max_step=0.03,t_eval=xs)
    u=xs*s.y[0]; Mm=np.column_stack([np.sin(2*np.pi*xs),np.cos(2*np.pi*xs)])
    (A,B),*_=np.linalg.lstsq(Mm,u,rcond=None); return B
Es=np.linspace(3,25,110); Bs=np.array([B_clean(E) for E in Es]); roots=[]
for i in range(len(Es)-1):
    if Bs[i]*Bs[i+1]<0:
        a,b=Es[i],Es[i+1]
        for _ in range(22):
            m=(a+b)/2
            if B_clean(a)*B_clean(m)<0: b=m
            else: a=m
        roots.append((a+b)/2)
roots=np.array(roots)
print(f"§463c clean(X=45). #roots[3,25]={len(roots)}, #zeros<25={int(np.sum(Z<25))}")
for k in range(min(12,len(roots))):
    gz=Z[k] if k<len(Z) else float('nan')
    print(f"  E={roots[k]:8.4f}  gamma={gz:8.4f}  ratio={roots[k]/gz if gz==gz else float('nan'):.4f}")

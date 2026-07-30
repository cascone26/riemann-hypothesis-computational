#!/usr/bin/env python3
"""§464b: sharpen — X=40, finer step, extract B over 8 periods. Confirm the match TIGHTENS (rules out coincidence)."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy')
lam=np.sqrt(2.0); FP=(2*np.pi*lam)**2; w=2*np.pi*lam; sp=lam
def Bval(gamma,X=40.0):
    alpha=-(gamma/2.0)**2; dxi0=-(FP*2-alpha)/(2*sp)
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-alpha)*y[0])/(2-x*x)]
    xs=np.linspace(X-6,X,900)
    s=solve_ivp(rhs,[sp+1e-5,X],[1.0,dxi0],rtol=1e-10,atol=1e-12,max_step=0.006,t_eval=xs)
    u=xs*s.y[0]; M=np.column_stack([np.sin(w*xs),np.cos(w*xs)])
    (A,B),*_=np.linalg.lstsq(M,u,rcond=None); return B
gs=np.linspace(12,60,300); Bs=np.array([Bval(g) for g in gs]); roots=[]
for i in range(len(gs)-1):
    if Bs[i]*Bs[i+1]<0:
        a,b=gs[i],gs[i+1]
        for _ in range(30):
            m=(a+b)/2
            if Bval(a)*Bval(m)<0: b=m
            else: a=m
        roots.append((a+b)/2)
roots=np.array(roots); zin=Z[(Z>=12)&(Z<=60)]
print("§464b sharpened (X=40). operator-γ vs zeta:")
print(f"{'k':>3} {'γ(op)':>10} {'γ(zeta)':>10} {'diff':>8}")
for k in range(min(len(roots),len(zin),16)):
    print(f"{k+1:>3} {roots[k]:10.4f} {zin[k]:10.4f} {roots[k]-zin[k]:+8.4f}")
print(f"#roots={len(roots)}, #zeros in [12,60]={len(zin)}; mean|diff|={np.mean(np.abs(roots[:min(len(roots),len(zin))]-zin[:min(len(roots),len(zin))])):.4f}")

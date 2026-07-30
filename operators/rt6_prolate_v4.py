#!/usr/bin/env python3
"""§464c: high-accuracy integration (DOP853, tight tol) over MODERATE X to beat phase-accumulation error.
Test if the operator-spectrum → zeta zeros tightens toward <1% (the paper's Fig-2 precision)."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy')
lam=np.sqrt(2.0); FP=(2*np.pi*lam)**2; w=2*np.pi*lam; sp=lam
def Bval(gamma,X):
    alpha=-(gamma/2.0)**2; dxi0=-(FP*2-alpha)/(2*sp)
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-alpha)*y[0])/(2-x*x)]
    xs=np.linspace(X-3/lam,X,700)   # last ~3 periods (period=1/√2)
    s=solve_ivp(rhs,[sp+1e-6,X],[1.0,dxi0],method='DOP853',rtol=1e-12,atol=1e-13,t_eval=xs)
    u=xs*s.y[0]; M=np.column_stack([np.sin(w*xs),np.cos(w*xs)])
    (A,B),*_=np.linalg.lstsq(M,u,rcond=None); return B
for X in [12.0, 16.0]:
    gs=np.linspace(12,40,240); Bs=np.array([Bval(g,X) for g in gs]); roots=[]
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1]<0:
            a,b=gs[i],gs[i+1]
            for _ in range(34):
                m=(a+b)/2
                if Bval(a,X)*Bval(m,X)<0: b=m
                else: a=m
            roots.append((a+b)/2)
    roots=np.array(roots); zin=Z[(Z>=12)&(Z<=40)]
    nn=min(len(roots),len(zin))
    md=np.mean(np.abs(roots[:nn]-zin[:nn])); mr=np.mean(np.abs(roots[:nn]-zin[:nn])/zin[:nn])*100
    print(f"X={X}: {len(roots)} roots vs {len(zin)} zeros; mean|diff|={md:.4f}, mean%err={mr:.2f}%")
    for k in range(min(8,nn)):
        print(f"   γ_op={roots[k]:8.4f}  γ_zeta={zin[k]:8.4f}  %={100*(roots[k]-zin[k])/zin[k]:+.2f}")

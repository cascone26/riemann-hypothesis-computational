#!/usr/bin/env python3
"""§464d: convergence test — DOP853, X ∈ {25,40,60}. If mean%err → 0 as X grows, the operator spectrum
converges to the zeta zeros (exact realization). Track the first ~8 zeros' convergence + Richardson extrapolate."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy')
lam=np.sqrt(2.0); FP=(2*np.pi*lam)**2; w=2*np.pi*lam; sp=lam
def Bval(gamma,X):
    alpha=-(gamma/2.0)**2; dxi0=-(FP*2-alpha)/(2*sp)
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-alpha)*y[0])/(2-x*x)]
    xs=np.linspace(X-3/lam,X,600)
    s=solve_ivp(rhs,[sp+1e-6,X],[1.0,dxi0],method='DOP853',rtol=1e-12,atol=1e-13,t_eval=xs)
    u=xs*s.y[0]; M=np.column_stack([np.sin(w*xs),np.cos(w*xs)])
    (A,B),*_=np.linalg.lstsq(M,u,rcond=None); return B
def roots_at(X):
    gs=np.linspace(12,35,190); Bs=np.array([Bval(g,X) for g in gs]); r=[]
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1]<0:
            a,b=gs[i],gs[i+1]
            for _ in range(32):
                m=(a+b)/2
                if Bval(a,X)*Bval(m,X)<0: b=m
                else: a=m
            r.append((a+b)/2)
    return np.array(r)
zin=Z[(Z>=12)&(Z<=35)]
res={}
for X in [25.0,40.0,60.0]:
    r=roots_at(X); nn=min(len(r),len(zin)); res[X]=r
    print(f"X={X:4.0f}: {len(r)} roots; mean%err(first {nn})={np.mean(np.abs(r[:nn]-zin[:nn])/zin[:nn])*100:.2f}%")
# Richardson: assume err ~ C/X; extrapolate roots to X->inf from X=40,60
r40,r60=res[40.0],res[60.0]; nn=min(len(r40),len(r60),len(zin))
rich=r60[:nn]+(r60[:nn]-r40[:nn])*(60)/(60-40)*0  # simple: if err~C/X, r(X)=r_inf+C/X => r_inf=(60*r60-40*r40)/20
rich=(60*r60[:nn]-40*r40[:nn])/20
print(f"\nRichardson extrapolation (X->inf, assuming err~C/X):")
print(f"{'k':>3} {'r40':>9} {'r60':>9} {'r_inf':>9} {'zeta':>9} {'%err':>7}")
for k in range(min(8,nn)):
    print(f"{k+1:>3} {r40[k]:9.4f} {r60[k]:9.4f} {rich[k]:9.4f} {zin[k]:9.4f} {100*(rich[k]-zin[k])/zin[k]:+7.2f}")
print(f"Richardson mean%err={np.mean(np.abs(rich[:min(8,nn)]-zin[:min(8,nn)])/zin[:min(8,nn)])*100:.2f}%")

#!/usr/bin/env python3
"""§463b (shooting on |x|>1, the correct region for negative eigenvalues). Regular solution at x=1⁺ (Frobenius):
ξ(1)=1, ξ'(1)=(μ-4π²)/2, φ(1)=0. Integrate ξ''=[2xξ'+((2π)²x²-μ)ξ]/(1-x²) on [1,X_∞]. At ∞ the solution ~
[A sin(2πx)+B cos(2πx)]/x; the W_sa condition (eq17) is B=0 (pure sin). Find μ<0 where B=0; check √|μ|=γ_k."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy'); FP=(2*np.pi)**2
def B_at_inf(E, X=40.0):
    mu=-E*E
    def rhs(x,y):
        return [y[1], (2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    e=1e-4
    y0=[1.0,(mu-FP)/2.0]
    s=solve_ivp(rhs,[1+e,X],y0,rtol=1e-9,atol=1e-11,max_step=0.02,dense_output=True)
    # extract A,B from u=x*ξ ≈ A sin(2πx)+B cos(2πx) near X: sample 2 points
    x1=X-0.5; x2=X-0.25
    u1=x1*s.sol(x1)[0]; u2=x2*s.sol(x2)[0]
    # solve [sin(2πx1) cos(2πx1); sin(2πx2) cos(2πx2)][A;B]=[u1;u2]
    M=np.array([[np.sin(2*np.pi*x1),np.cos(2*np.pi*x1)],[np.sin(2*np.pi*x2),np.cos(2*np.pi*x2)]])
    A,B=np.linalg.solve(M,[u1,u2])
    return B, A
Es=np.linspace(3,32,300)
Bs=np.array([B_at_inf(E)[0] for E in Es])
roots=[]
for i in range(len(Es)-1):
    if Bs[i]*Bs[i+1]<0:
        a,b=Es[i],Es[i+1]
        for _ in range(35):
            m=(a+b)/2
            if B_at_inf(a)[0]*B_at_inf(m)[0]<0: b=m
            else: a=m
        roots.append((a+b)/2)
roots=np.array(roots)
print(f"§463b shooting on |x|>1, B=0 (pure-sin) condition. E=√|μ| vs zeta zeros:")
print(f"{'k':>3} {'E=√|μ|':>10} {'γ_k':>10} {'ratio':>7}")
for k in range(min(20,len(roots))):
    gz=Z[k] if k<len(Z) else np.nan
    print(f"{k+1:>3} {roots[k]:10.4f} {gz:10.4f} {roots[k]/gz if gz==gz else float('nan'):7.4f}")
print(f"\n#roots in E∈[3,32]={len(roots)}, #zeros<32={int(np.sum(Z<32))}")
print("reading: E-values matching γ_k => CORRECT realization of the zeros by the CCM operator (built correctly!).")

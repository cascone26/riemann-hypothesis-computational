#!/usr/bin/env python3
"""§464 (CORRECTED per paper §6): λ=√2 (not 1!). W_{√2}=-∂_x[(2-x²)∂_x]+8π²x². Negative eigenvalues α_k<0.
The zeros are γ_k = 2√|α_k| (Dirac 2D: eigenvalues ±2√α, imag parts 2√|α| = the zeros, Fig 2). Asymptotic
frequency ω=2π√2. Regular at x=√2: ξ(√2)=1, ξ'(√2)=(α-16π²)/(2√2). Pure-sin B=0 at ∞. Shoot for B=0, parametrized
by γ (α=-(γ/2)²), and compare to the zeta zeros DIRECTLY."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy')
lam=np.sqrt(2.0); lam2=2.0; FP=(2*np.pi*lam)**2   # (2πλ)²=8π²
w=2*np.pi*lam                                       # asymptotic frequency 2π√2
sp=lam                                              # singular point √2
def Bval(gamma, X=22.0):
    alpha=-(gamma/2.0)**2
    dxi0=(alpha-FP*lam2)/(2*lam*np.sqrt(lam2))       # wait: regularity: 2*sp*ξ' + (FP*sp²-α)ξ=0 ; sp²=2
    dxi0=-(FP*lam2-alpha)/(2*sp)                      # ξ'(√2) = -(8π²·2-α)/(2√2)
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-alpha)*y[0])/(lam2-x*x)]
    xs=np.linspace(X-4,X,500)
    s=solve_ivp(rhs,[sp+1e-4,X],[1.0,dxi0],rtol=1e-9,atol=1e-11,max_step=0.01,t_eval=xs)
    u=xs*s.y[0]; M=np.column_stack([np.sin(w*xs),np.cos(w*xs)])
    (A,B),*_=np.linalg.lstsq(M,u,rcond=None); return B
# scan gamma over the zeta-zero range, find B=0
gs=np.linspace(10,45,240); Bs=np.array([Bval(g) for g in gs])
roots=[]
for i in range(len(gs)-1):
    if Bs[i]*Bs[i+1]<0:
        a,b=gs[i],gs[i+1]
        for _ in range(28):
            m=(a+b)/2
            if Bval(a)*Bval(m)<0: b=m
            else: a=m
        roots.append((a+b)/2)
roots=np.array(roots)
print(f"§464 CORRECTED (λ=√2, γ=2√|α|, freq 2π√2). Spectrum-γ vs zeta zeros:")
print(f"{'k':>3} {'γ (operator)':>13} {'γ_k (zeta)':>11} {'diff':>8} {'ratio':>7}")
zin=Z[(Z>=10)&(Z<=45)]
for k in range(min(len(roots),len(zin),22)):
    print(f"{k+1:>3} {roots[k]:13.4f} {zin[k]:11.4f} {roots[k]-zin[k]:+8.4f} {roots[k]/zin[k]:7.4f}")
print(f"\n#operator-roots in [10,45]={len(roots)}, #zeros in [10,45]={len(zin)}")
print("=> if γ(operator) ≈ γ_k(zeta), I've BUILT THE CCM OPERATOR CORRECTLY and it realizes the zeros.")

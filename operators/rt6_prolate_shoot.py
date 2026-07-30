#!/usr/bin/env python3
"""§463 (CORRECT extension via shooting): the W_sa eigenfunctions must be REGULAR at x=1 (flux φ=(1-x²)ξ'→0).
Even solution from x=0: ξ(0)=1, φ(0)=0. System: ξ'=φ/(1-x²), φ'=((2π)²x²-μ)ξ. For generic μ, φ(1)≠0 ⟹ ξ' blows
up at 1 (singular, NOT in domain). The special μ where φ(1)=0 are REGULAR (candidate eigenvalues). Compute
φ(1⁻;μ) vs μ<0, find its zeros, check if √|μ| = the zeta zeros γ_k. This tests the correct-extension realization."""
import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('zeros500.npy'); FP=(2*np.pi)**2
def flux_at_1(mu):
    # integrate (ξ,φ) from 0 to 1-eps; return φ(1-eps)
    def rhs(x,y): return [y[1]/(1-x*x), (FP*x*x-mu)*y[0]]
    s=solve_ivp(rhs,[0,1-1e-6],[1.0,0.0],rtol=1e-10,atol=1e-12,dense_output=False,max_step=0.001)
    return s.y[1,-1]   # φ near 1
# scan mu<0, find sign changes of φ(1;μ)
mus=-np.linspace(1,700,7000)**2/ (2*np.pi)**0  # will refine; use mu = -E^2 with E in [1,~26]
Es=np.linspace(3,60,4000); mus=-Es**2
fp=np.array([flux_at_1(m) for m in mus])
# find sign changes -> roots
roots=[]
for i in range(len(mus)-1):
    if fp[i]*fp[i+1]<0:
        # linear interp for E
        E0,E1=Es[i],Es[i+1]; f0,f1=fp[i],fp[i+1]
        Er=E0-f0*(E1-E0)/(f1-f0); roots.append(Er)
roots=np.array(roots)
print(f"§463 shooting: μ where even-from-0 solution is REGULAR at x=1 (φ(1)=0). √|μ|=E values (candidate zeros):")
print(f"{'k':>3} {'E=√|μ| (regular)':>16} {'γ_k (zeta)':>10} {'ratio':>7}")
for k in range(min(15,len(roots))):
    gz=Z[k] if k<len(Z) else np.nan
    print(f"{k+1:>3} {roots[k]:16.4f} {gz:10.4f} {roots[k]/gz if gz==gz else np.nan:7.4f}")
print(f"\n#roots found in E∈[3,60]: {len(roots)}. Compare to #zeros<60 = {int(np.sum(Z<60))}.")
print("""
reading: if the regular-at-1 E-values match the zeta zeros γ_k, the correct-extension shooting realizes the zeros
— I've built the CCM operator correctly. If they don't match but have the right COUNT/density, the regularity-at-1
condition is necessary but the ∞ (B=0) condition also matters (need both). If count is way off, the shooting setup
needs the discontinuity-allowed matching (full W_sa).
""")

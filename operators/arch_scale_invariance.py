#!/usr/bin/env python3
"""Cycle 9: verify the analytic decomposition  lambda_arch(a) = mu1 - log(a)  (arch-only, no primes).
Suzuki (2.4): La(v) = 1/4 ∫∫|v(x)-v(y)|^2/|x-y| - 1/2 ∫ log(a^2-x^2)|v|^2, over (-a,a).
Split log(a^2-x^2)=2log a + log(1-(x/a)^2). H^{1/2} Dirichlet energy is SCALE-INVARIANT in 1D, and the
edge term 1/2∫log(1-(x/a)^2)|v|^2 is scale-invariant under x=a*u. So La(v)=[scale-inv arch]-log a·||v||^2
=> lambda_arch(a)=mu1 - log a, mu1 a-INDEPENDENT. PREDICTS: lambda_arch(a)+log a = mu1 = const ~0.2368.
Verify-before-ship: compute lambda_arch at several a directly on (-a,a), check +log a is flat = mu1."""
import numpy as np
from scipy.linalg import eigh
def lam_arch(a,M=4000,K=100):
    x=np.linspace(-a,a,M+1)[1:-1]; h=2*a/M
    kk=np.arange(1,K+1); Vk=np.sin(np.pi*np.outer((x+a),kk)/(2*a)); w=np.full(len(x),h)
    xr=x[:,None]-x[None,:]; Kk=np.where(np.abs(xr)>1e-9,1.0/np.abs(xr),0.0); Krow=(Kk*w[None,:]).sum(1)
    # 1/4 ∫∫|v(x)-v(y)|^2/|x-y| = 1/4[ 2∫|v|^2 Krow - 2∫∫ v(x)v(y)/|x-y| ]  (quadratic form on basis)
    jump=0.25*((Vk.T*(2*w*Krow))@Vk - 2*Vk.T@(w[:,None]*Kk*w[None,:])@Vk)
    kill=-0.5*(Vk.T*(w*np.log(a*a-x*x)))@Vk       # -1/2 ∫ log(a^2-x^2)|v|^2
    L=(jump+kill);L=(L+L.T)/2; Mass=Vk.T@(w[:,None]*Vk);Mass=(Mass+Mass.T)/2
    return eigh(L,Mass,eigvals_only=True)[0]
print(f"{'a':>5} {'lam_arch(a)':>12} {'+log a':>10} {'pred mu1-log a':>15}")
mu1=0.23676
for a in [0.7,1.0,1.5,2.0,3.0,4.0,6.0]:
    la=lam_arch(a); flat=la+np.log(a); pred=mu1-np.log(a)
    print(f"{a:5.1f} {la:+12.5f} {flat:+10.5f} {pred:+15.5f}")
print(f"\nIf '+log a' column is FLAT ~ {mu1} across all a => VERIFIED: lambda_arch(a)=mu1-log a analytically.")
print("Then FULL lambda_a = mu1 - log a + prime_form(a); RH <=> prime_form(a) >= log a - mu1 forall a.")

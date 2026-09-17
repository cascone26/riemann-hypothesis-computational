#!/usr/bin/env python3
"""Cycle 23: characterize the universal profile Phi(u)=lim a*M_full(u/a;a). Compute finely at a=4 (converged),
fit small-u law, locate well+zeros, and test whether Phi is the explicit-formula fluctuation (arch log - primes)."""
import numpy as np, suzuki_screw_lambda_a as S
from scipy.special import digamma
a=4.0
Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
t=np.linspace(-2*a,2*a,20001); g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])
hann=0.5*(1+np.cos(np.pi*t/(2*a))); gh=g*hann
u=np.linspace(0.05,10,2000); z=u/a
ghat=np.array([np.trapezoid(gh*np.cos(zz*t),t) for zz in z]); Phi=z*z*ghat*a
# small-u law: fit Phi ~ c*u^p for u<0.8
m=u<0.8; p=np.polyfit(np.log(u[m]),np.log(-Phi[m]),1)
print(f"small-u: Phi(u) ~ {-np.exp(p[1]):.4f} * u^{p[0]:.3f}  (u<0.8)")
# well + zeros
iwell=np.argmin(Phi); print(f"well: Phi_min={Phi[iwell]:.5f} at u={u[iwell]:.4f}")
sign=np.sign(Phi); zc=u[:-1][sign[:-1]!=sign[1:]]
print(f"zero-crossings (u): {np.round(zc[:8],3)}")
print(f"zero spacings: {np.round(np.diff(zc[:8]),3)}")
# compare Phi to candidate: arch symbol log(u)+C0 shifted vs the actual (is Phi the fluctuation only?)
C0=0.5772156649
print(f"\nPhi values: "+" ".join(f"{uu}:{np.interp(uu,u,Phi):+.4f}" for uu in [0.5,1,1.84,3,3.8,5,7,9]))
# does Phi integrate to 0? (mean-zero-ish / explicit-formula sum rule)
print(f"int Phi(u) du [0,10] = {np.trapezoid(Phi,u):+.4f}   int Phi/u^2 du = {np.trapezoid(Phi/u**2,u):+.4f}")

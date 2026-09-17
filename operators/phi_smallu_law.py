#!/usr/bin/env python3
"""Cycle 26: the TRUE small-u law of Phi_central(u)=lim a*M_full(u/a;a) (hard cutoff=exact). Governs large-a.
Fit Phi ~ c*u^p for small u; identify c (vs mu1=0.2368, log2/sqrt2=0.490, etc). Converged a=4."""
import numpy as np, suzuki_screw_lambda_a as S
a=4.0; Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
t=np.linspace(-2*a,2*a,60001); g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])
u=np.linspace(0.05,1.4,60); z=u/a
ghat=np.array([np.trapezoid(g*np.cos(zz*t),t) for zz in z]); Phi=z*z*ghat*a
# fit on several small-u windows to see the exponent stabilize
for hi in [0.4,0.6,0.9]:
    m=(u>0.1)&(u<hi); p=np.polyfit(np.log(u[m]),np.log(-Phi[m]),1)
    print(f"  u in (0.1,{hi}): Phi ~ {-np.exp(p[1]):+.4f} * u^{p[0]:.3f}")
# also fit Phi = c1*u + c2*u^2 (mixed)
m=u<1.0; c=np.polyfit(u[m],Phi[m],2)
print(f"  quad fit Phi ~ {c[0]:+.4f} u^2 {c[1]:+.4f} u {c[2]:+.4f}")
print(f"  samples: "+" ".join(f"u{uu}:{np.interp(uu,u,Phi):+.4f}" for uu in [0.1,0.2,0.4,0.7,1.0,1.3]))
# candidate constants
import math
print(f"\n  mu1=0.2368  log2/sqrt2={math.log(2)/math.sqrt(2):.4f}  1/(2pi)={1/(2*math.pi):.4f}  A=0.7075")

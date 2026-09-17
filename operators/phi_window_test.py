#!/usr/bin/env python3
"""Cycle 23 verify: is Phi(u) real screw content or a Hann-window artifact? Recompute a*M_full(u/a;a) with
THREE windows (Hann, Gaussian exp(-(t/a)^2*k), hard/none) at a=4. If well+zeros PERSIST => real; if they
move with the window => Phi is the window's transform (artifact), deflating cycles 22-23."""
import numpy as np, suzuki_screw_lambda_a as S
a=4.0; Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
t=np.linspace(-2*a,2*a,20001); g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])
u=np.linspace(0.05,10,2000); z=u/a
wins={
 'hann':0.5*(1+np.cos(np.pi*t/(2*a))),
 'gauss':np.exp(-(t/a)**2*0.8),
 'tukey0.5':np.where(np.abs(t)<a, 1.0, 0.5*(1+np.cos(np.pi*(np.abs(t)-a)/a))),
}
print(f"{'window':>10} {'well u':>8} {'well val':>9} {'zeros u (first 4)':>26}")
for nm,win in wins.items():
    gh=g*win; ghat=np.array([np.trapezoid(gh*np.cos(zz*t),t) for zz in z]); Phi=z*z*ghat*a
    iw=np.argmin(Phi); sign=np.sign(Phi); zc=u[:-1][sign[:-1]!=sign[1:]]
    print(f"{nm:>10} {u[iw]:8.3f} {Phi[iw]:9.4f}   {np.round(zc[:4],3)}")
print("\nIf well u & zeros are the SAME across windows => real screw content. If they MOVE => window artifact.")

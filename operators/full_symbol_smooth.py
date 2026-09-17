#!/usr/bin/env python3
"""Cycle 20: M_full(z;a) with a SMOOTH (Hann) window to kill the |t|=2a Gibbs artifact; profile at moderate z
(where the ground-state mass lives, z-mean~0.5-1.3). Determine the REAL sign structure of the complete symbol."""
import numpy as np, suzuki_screw_lambda_a as S
def prof(a):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    t=np.linspace(-2*a,2*a,8001); dt=t[1]-t[0]
    g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])
    hann=0.5*(1+np.cos(np.pi*t/(2*a)))   # smooth taper to 0 at |t|=2a
    gh=g*hann
    z=np.linspace(0.05,12,240)
    ghat=np.array([np.trapezoid(gh*np.cos(zz*t),t) for zz in z])
    return z,z*z*ghat
print("M_full(z;a) [Hann-windowed] at moderate z — sign structure where the mass lives:")
for a in [1.0,1.5,2.0,2.5]:
    z,M=prof(a)
    neg=z[M<0]
    firstneg = neg[0] if len(neg) else None
    # sample values
    samp=[(zz,np.interp(zz,z,M)) for zz in [0.5,1,2,3,5,8]]
    print(f" a={a}: min={M.min():+.3f}@z={z[np.argmin(M)]:.2f}  M<0 on z in "
          f"[{firstneg if firstneg else '--'}..]  samples "+" ".join(f"z{zz:.0f}:{v:+.2f}" for zz,v in samp))
print("\nIf M_full>=0 for small-moderate z (where mass sits) => positivity plausibly lift-free; if M_full<0")
print("right where the mass concentrates => the negativity is real & the lift/uncertainty is essential.")

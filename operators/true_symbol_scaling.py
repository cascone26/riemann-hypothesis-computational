#!/usr/bin/env python3
"""Cycle 24: does the TRUE (hard-cutoff = exact kernel support) symbol collapse under u=za? The kernel g(x-y)
is only evaluated at |t|<2a on (-a,a), so int_{-2a}^{2a} is EXACT (no window). Test a*M_full(u/a;a) collapse
with the HARD cutoff (fine grids), and compare the true-symbol shape to the Hann-contaminated cycle-22 Phi."""
import numpy as np, suzuki_screw_lambda_a as S
def prof(a,uu,NT=40001):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    t=np.linspace(-2*a,2*a,NT); g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])   # HARD, exact support
    z=uu/a; ghat=np.array([np.trapezoid(g*np.cos(zz*t),t) for zz in z])
    return z*z*ghat*a
uu=np.array([0.5,1.0,1.5,1.84,2.5,3.14,4.0,6.0,10.0,20.0,35.0])
print("TRUE hard-cutoff symbol a*M_full(u/a;a):")
print(f"{'u':>6} "+" ".join(f"{'a='+str(a):>10}" for a in [1.0,2.0,3.0]))
cols={a:prof(a,uu) for a in [1.0,2.0,3.0]}
for i,u in enumerate(uu):
    print(f"{u:6.2f} "+" ".join(f"{cols[a][i]:+10.4f}" for a in [1.0,2.0,3.0]))
M=np.array([cols[a] for a in cols]); rel=np.nanstd(M,axis=0)/(np.abs(np.nanmean(M,0))+1e-6)
print(f"\nrel-spread across a per u: {np.round(rel,3)}")
print("small-u (u<=1.84): collapse? tail (u>=4): collapse? If tail does NOT collapse => cycle-22 universal")
print("profile was Hann-driven; only the small-u/central part is real.")

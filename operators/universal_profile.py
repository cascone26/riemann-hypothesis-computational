#!/usr/bin/env python3
"""Cycle 22: does M_full(u/a;a)*a collapse to a universal profile Phi(u), u=z*a? (dip depth*a~0.117 const,
dip loc z*a~1.76 const suggest yes.) Compute M_full*a vs u for several a; if they overlay => a-indep inequality."""
import numpy as np, suzuki_screw_lambda_a as S
def prof(a,uu):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    t=np.linspace(-2*a,2*a,8001); g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])
    hann=0.5*(1+np.cos(np.pi*t/(2*a))); gh=g*hann
    z=uu/a
    ghat=np.array([np.trapezoid(gh*np.cos(zz*t),t) for zz in z])
    return z*z*ghat*a   # M_full * a
uu=np.linspace(0.2,6,40)
print("u=z*a  |  M_full(u/a;a)*a  for a=1,1.5,2,2.5,3  (collapse => universal Phi(u)):")
cols={a:prof(a,uu) for a in [1.0,1.5,2.0,2.5,3.0]}
print(f"{'u':>5} "+" ".join(f"{'a='+str(a):>9}" for a in cols))
for i in range(0,len(uu),3):
    print(f"{uu[i]:5.2f} "+" ".join(f"{cols[a][i]:+9.4f}" for a in cols))
# quantify collapse: std across a at each u, relative to mean
M=np.array([cols[a] for a in cols])
rel=np.nanstd(M,axis=0)/(np.abs(np.nanmean(M,axis=0))+1e-9)
print(f"\nmean rel-spread across a (u in [1,4]): {np.nanmean(rel[(uu>=1)&(uu<=4)]):.3f}  (small => good collapse)")
print("dip of the a=3 profile: min=%.4f at u=%.2f"%(cols[3.0].min(),uu[np.argmin(cols[3.0])]))

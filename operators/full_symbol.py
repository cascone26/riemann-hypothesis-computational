#!/usr/bin/env python3
"""Cycle 19: the COMPLETE localized Weil multiplier. Form = int int g(x-y) w'(x)w'(y) = (1/2pi)int ghat(z) z^2 |what|^2.
So M_full(z;a) = z^2 * ghat_windowed(z),  ghat_windowed(z)=int_{-2a}^{2a} g(t) cos(zt) dt  (g even, |x-y|<2a).
Includes EVERYTHING (arch+pole+prime+r'') via the full screw g. Check <M_full>_dens = lam_full~0, then sign."""
import numpy as np, suzuki_screw_lambda_a as S, suzuki_exact_Q as EQ
def run(a,M=3000,K=50,nz=3000):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    # full screw g on t in (-2a,2a)
    t=np.linspace(-2*a,2*a,8001); dt=t[1]-t[0]
    g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])
    z=np.linspace(0.02,40,nz)
    ghat=np.array([np.trapezoid(g*np.cos(zz*t),t) for zz in z])
    Mfull=z*z*ghat
    # ground state density in Fourier (same as cycle 18)
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h; idx=np.arange(M+1)
    Gf=g[0]*0  # placeholder
    gg=np.array([S.screw_g(d,Lam,Nmax) for d in diffs]); Gf=gg[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a); Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2
    ev,V=np.linalg.eigh(Qf); c0=V[:,0]; lam=ev[0]/a
    Sh=np.array([EQ.shat(zz,a,kk) for zz in z]); dens=np.abs(Sh@c0)**2
    Mavg=np.trapezoid(Mfull*dens,z)/np.trapezoid(dens,z)/a   # density-avg symbol (/a for lam normalization)
    negmass=np.trapezoid(dens*(Mfull<0),z)/np.trapezoid(dens,z)
    minM=Mfull.min(); zmin=z[np.argmin(Mfull)]
    return lam,Mavg,negmass,minM,zmin
print(f"{'a':>5} {'lam_full':>9} {'<M_full>/a':>10} {'negmass':>8} {'min M_full':>11} {'z@min':>7}")
for a in [1.0,1.5,2.0,2.5]:
    lam,Mavg,nm,minM,zmin=run(a)
    print(f"{a:5.1f} {lam:+9.4f} {Mavg:+10.4f} {nm:8.3f} {minM:+11.3f} {zmin:7.2f}")
print("\n<M_full>/a should ~ lam_full (validates the complete symbol). Then min M_full & where: is the FULL")
print("symbol >=0 away from z=0 (=> lift-free positivity via 1759) or still negative (genuine crux)?")

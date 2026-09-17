#!/usr/bin/env python3
"""Cycle 21 redo: HARD-window symbol (quantitatively = operator, validated <M>/a=lam_full~0) but integrate
only z<Zcut to exclude the z~35 Gibbs artifact. Check pos/neg parts of the Rayleigh integral and their scaling.
Also report the full-range integral to confirm it ~ lam_full~0 (so the artifact contributes ~0 net)."""
import numpy as np, suzuki_screw_lambda_a as S, suzuki_exact_Q as EQ
def run(a,M=3000,K=50,Zcut=15):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    t=np.linspace(-2*a,2*a,8001); g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])   # HARD window
    z=np.linspace(0.02,45,9000)
    ghat=np.array([np.trapezoid(g*np.cos(zz*t),t) for zz in z]); Mf=z*z*ghat
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h; idx=np.arange(M+1)
    gg=np.array([S.screw_g(d,Lam,Nmax) for d in diffs]); Gf=gg[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a); Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2; ev,V=np.linalg.eigh(Qf); c0=V[:,0]
    Sh=np.array([EQ.shat(zz,a,kk) for zz in z]); dens=np.abs(Sh@c0)**2; norm=np.trapezoid(dens,z)
    full=np.trapezoid(Mf*dens,z)/norm/a
    m=z<Zcut
    Ip=np.trapezoid(np.where(Mf[m]>0,Mf[m],0)*dens[m],z[m])/norm/a
    Im=np.trapezoid(np.where(Mf[m]<0,Mf[m],0)*dens[m],z[m])/norm/a
    return full,Ip,Im
print(f"{'a':>5} {'full/a~lam':>11} {'I+(z<15)':>9} {'I-(z<15)':>9} {'|I-|*a':>8} {'|I-|*a^2':>9}")
for a in [1.0,1.5,2.0,2.5,3.0]:
    full,Ip,Im=run(a)
    print(f"{a:5.1f} {full:+11.5f} {Ip:+9.4f} {Im:+9.4f} {abs(Im)*a:8.4f} {abs(Im)*a*a:9.4f}")
print("\nfull/a should ~0 (validates). I+/I- over z<15 (artifact-free): scaling of |I-| tests the mechanism.")

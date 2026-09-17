#!/usr/bin/env python3
"""Cycle 21: is the lambda_full~0 cancellation PROTECTED (both pos/neg parts shrink) or DELICATE (both large)?
Decompose the Fourier Rayleigh integral of the COMPLETE symbol on the actual ground state:
  I+ = int_{M>0} M_full |what|^2,  I- = int_{M<0} M_full |what|^2,  lam_full ~ (I+ + I-)/(||.||*a).
Track |I-| vs a. Mechanism (dip retreats into mean-zero zone) predicts |I-| shrinks as ~1/a^2."""
import numpy as np, suzuki_screw_lambda_a as S, suzuki_exact_Q as EQ
def run(a,M=3000,K=50):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    t=np.linspace(-2*a,2*a,8001); g=np.array([S.screw_g(tt,Lam,Nmax) for tt in t])
    hann=0.5*(1+np.cos(np.pi*t/(2*a))); gh=g*hann
    z=np.linspace(0.02,20,4000)
    ghat=np.array([np.trapezoid(gh*np.cos(zz*t),t) for zz in z]); Mf=z*z*ghat
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h; idx=np.arange(M+1)
    gg=np.array([S.screw_g(d,Lam,Nmax) for d in diffs]); Gf=gg[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a); Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2; ev,V=np.linalg.eigh(Qf); c0=V[:,0]
    Sh=np.array([EQ.shat(zz,a,kk) for zz in z]); dens=np.abs(Sh@c0)**2
    norm=np.trapezoid(dens,z)
    Ip=np.trapezoid(np.where(Mf>0,Mf,0)*dens,z)/norm
    Im=np.trapezoid(np.where(Mf<0,Mf,0)*dens,z)/norm
    return Ip,Im,ev[0]/a
print(f"{'a':>5} {'I+ (pos)':>9} {'I- (neg)':>9} {'|I-|':>8} {'|I-|*a':>8} {'|I-|*a^2':>9} {'lam~(I++I-)':>11}")
for a in [1.0,1.5,2.0,2.5,3.0]:
    Ip,Im,lam=run(a)
    print(f"{a:5.1f} {Ip:+9.4f} {Im:+9.4f} {abs(Im):8.4f} {abs(Im)*a:8.4f} {abs(Im)*a*a:9.4f} {Ip+Im:+11.4f}")
print("\nIf |I-|*a^2 ~ const => neg part shrinks as 1/a^2 (mechanism holds, positivity increasingly protected).")
print("If |I-| ~ const or grows => delicate cancellation persists (mechanism weaker than hoped).")

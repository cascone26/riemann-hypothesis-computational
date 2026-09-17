#!/usr/bin/env python3
"""Cycle 18: where does the actual lambda_a ground state live in Fourier vs where M(z;a) is negative?
Full ground eigvec c0 (screw frame) -> w-hat(z)=sum c0_k shat(z,a,k). Mass |w-hat|^2 distribution vs M(z;a).
The z-band carrying the mass = where the uncertainty inequality is tight (where a sharp majorant must be built)."""
import numpy as np, suzuki_screw_lambda_a as S, suzuki_exact_Q as EQ
C0=0.5772156649; A2p1=np.log(2*np.pi)+C0
def run(a,M=3000,K=50):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h; idx=np.arange(M+1)
    gf=np.array([S.screw_g(d,Lam,Nmax) for d in diffs]); Gf=gf[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a); Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2
    ev,V=np.linalg.eigh(Qf); c0=V[:,0]; lam=ev[0]/a
    # w-hat(z) = sum_k c0_k * shat(z,a,k)  (shat uses kk index)
    z=np.linspace(0.01,50,4000)
    Sh=np.array([EQ.shat(zz,a,kk) for zz in z])   # (nz,K)
    what=Sh@c0
    dens=np.abs(what)**2; dens/=np.trapezoid(dens,z)
    # M(z;a)
    ps=np.zeros_like(z)
    for n in range(2,min(Nmax,int(np.exp(2*a)))+1):
        if Lam[n]>0 and np.log(n)<2*a: ps+=Lam[n]/np.sqrt(n)*np.cos((np.log(n)/a)*z)
    Msym=np.log(np.abs(z))+C0-np.log(a)-A2p1-2*ps
    # mass-weighted M (should ~ lam*something) and where mass sits
    zmean=np.trapezoid(z*dens,z); zpeak=z[np.argmax(dens)]
    Mneg_frac=np.trapezoid(dens*(Msym<0),z)
    return lam,zmean,zpeak,Mneg_frac,np.trapezoid(dens*Msym,z)
print(f"{'a':>5} {'lam_full':>9} {'z-mean':>7} {'z-peak':>7} {'mass frac where M<0':>19} {'<M>_dens':>9}")
for a in [1.0,1.5,2.0,2.5]:
    lam,zm,zp,mf,Mw=run(a)
    print(f"{a:5.1f} {lam:+9.4f} {zm:7.2f} {zp:7.2f} {mf:19.3f} {Mw:+9.4f}")
print("\nz-band of the mass = where the inequality is tight. mass-frac where M<0: how much sits in the")
print("negative-symbol region (the lift must overcome). <M>_dens = density-averaged multiplier.")

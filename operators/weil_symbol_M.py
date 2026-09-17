#!/usr/bin/env python3
"""Cycle 14: the localized Weil Fourier multiplier M(z;a) from Suzuki (4.5)+(4.6)+Wiener-Khinchin.
M(z;a)=log|z|+C0 -log a -(2A+1) -2*Sum_{n<=e^2a} Lambda(n)/sqrt(n) cos((log n/a) z)   [drop r'' for leading symbol M0]
Since w is supported in (-1,1), lambda_a >= ess-inf_z M(z;a) (RIGOROUS gap-free lower bound; support/uncertainty
lifts lambda_a above the symbol min). Compute M0 min vs a; compare to screw lambda_full(~0) and lambda_arch."""
import numpy as np, suzuki_screw_lambda_a as S
C0=0.5772156649; A2p1=np.log(2*np.pi)+C0   # 2A+1
def Msym(a,zmax=60,nz=6000):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,10_000_000))
    z=np.linspace(1e-3,zmax,nz)
    base=np.log(np.abs(z))+C0-np.log(a)-A2p1
    ps=np.zeros_like(z)
    for n in range(2,min(Nmax,int(np.exp(2*a)))+1):
        if Lam[n]>0 and np.log(n)<2*a:
            ps+=Lam[n]/np.sqrt(n)*np.cos((np.log(n)/a)*z)
    M=base-2*ps
    return z,M
print(f"{'a':>5} {'min M0(z)':>10} {'argmin z':>9} {'M0(z=0+)':>10} {'M0(zbig)':>10}")
larch={0.4:-0.072,0.7:-0.727,1.0:-1.312,1.5:-2.621,2.0:-4.897,2.5:-8.925,3.0:-15.907}
for a in [0.4,0.7,1.0,1.5,2.0,2.5,3.0]:
    z,M=Msym(a); i=np.argmin(M)
    print(f"{a:5.1f} {M[i]:+10.3f} {z[i]:9.2f} {M[0]:+10.3f} {M[-1]:+10.3f}   (screw lam_arch={larch[a]:+.2f}, lam_full~0)")
print("\nlambda_a >= min_z M0 (minus r'' correction). If min M0 >> lambda_arch, the support-lift is large;")
print("track whether min M0 is closer to lambda_full(~0) [=> symbol nearly governs] or to lambda_arch.")

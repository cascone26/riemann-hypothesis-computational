#!/usr/bin/env python3
"""Cycle 25: decompose the TRUE universal central symbol (u<=2.5, hard cutoff) into ARCH-only vs PRIME-only.
a*M(u/a;a) with g_arch (screw_g_no_primes) and g_prime = g - g_arch, separately. Locates the shallow-dip origin
and gets toward a closed form (arch = FT of the scale-invariant screw core)."""
import numpy as np, suzuki_screw_lambda_a as S
def prof(a,uu,gfun,NT=40001):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    t=np.linspace(-2*a,2*a,NT); g=np.array([gfun(tt,Lam,Nmax) for tt in t])
    z=uu/a; ghat=np.array([np.trapezoid(g*np.cos(zz*t),t) for zz in z])
    return z*z*ghat*a
uu=np.array([0.5,1.0,1.5,1.84,2.5])
garch=lambda tt,L,N: S.screw_g_no_primes(tt,L,N)
gprime=lambda tt,L,N: S.screw_g(tt,L,N)-S.screw_g_no_primes(tt,L,N)
gfull=lambda tt,L,N: S.screw_g(tt,L,N)
print("central symbol decomposition (a=2, hard cutoff):")
print(f"{'u':>6} {'ARCH':>9} {'PRIME':>9} {'FULL':>9} {'arch+prime':>10}")
A=prof(2.0,uu,garch); P=prof(2.0,uu,gprime); F=prof(2.0,uu,gfull)
for i,u in enumerate(uu):
    print(f"{u:6.2f} {A[i]:+9.4f} {P[i]:+9.4f} {F[i]:+9.4f} {A[i]+P[i]:+10.4f}")
# universality of arch vs prime parts separately
A1=prof(1.0,uu,garch);A3=prof(3.0,uu,garch)
P1=prof(1.0,uu,gprime);P3=prof(3.0,uu,gprime)
print("\narch collapse (a=1,2,3) at u=1:  %.4f %.4f %.4f"%(A1[1],A[1],A3[1]))
print("prime collapse (a=1,2,3) at u=1: %.4f %.4f %.4f"%(P1[1],P[1],P3[1]))
print("\nIf arch is smooth+universal and prime carries the dip => the dip is the first-prime (n=2,log2) feature.")

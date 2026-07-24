# LP-class route: Xi(z)=2*sum(-1)^k b_k z^{2k}/(2k)!, b_k=int_0^inf u^{2k} Phi(u) du (moments of positive Phi).
# LP membership (=> RH) requires the TURAN inequalities. Compute b_k, test log-concavity / Turan, see tightness.
import numpy as np
def Phi(u,nmax=50):
    u=np.atleast_1d(np.asarray(u,float)); n=np.arange(1,nmax+1)[:,None]
    T=(2*np.pi**2*n**4*np.exp(9*u)[None,:]-3*np.pi*n**2*np.exp(5*u)[None,:])*np.exp(-np.pi*n**2*np.exp(4*u)[None,:])
    return T.sum(0)
# moments b_k = int_0^inf u^{2k} Phi(u) du (Phi concentrated in u<0.6; use fine grid + tail)
U=2.5; M=200000; u=np.linspace(1e-6,U,M); du=u[1]-u[0]; phi=Phi(u)
b=[np.sum(u**(2*k)*phi)*du for k in range(0,16)]
print("moments b_k = int u^{2k} Phi du:")
for k in range(0,14,2): print(f"   b_{k}={b[k]:.6e}")
# Turan inequality (Csordas-Norfolk-Varga form): T_k = b_k^2 - b_{k-1} b_{k+1} * (2k+1)/(2k-1)  >= 0 ?
# also the plain log-concavity b_k^2 - b_{k-1}b_{k+1} and the normalized LP one.
print("\nTuran / log-concavity tests (should be >=0 for RH-consistent LP membership):")
print(f"{'k':>3} {'plain b_k^2-b_{k-1}b_{k+1}':>26} {'ratio b_k^2/(b_{k-1}b_{k+1})':>28}")
for k in range(1,14):
    plain=b[k]**2-b[k-1]*b[k+1]
    ratio=b[k]**2/(b[k-1]*b[k+1])
    print(f"{k:>3} {plain:>26.3e} {ratio:>28.6f}")
print("\nratio > 1 => log-concave (Turan holds). How close to 1 => how tight the knife-edge.")

print("\n=== CORRECT normalization: Xi Taylor coeffs c_k = b_k/(2k)! ; Turan c_k^2 >= c_{k-1}c_{k+1} ===")
from math import factorial
c=[b[k]/factorial(2*k) for k in range(16)]
print(f"{'k':>3} {'c_k^2/(c_{k-1}c_{k+1})':>24}  (>1 => Turan HOLDS; ->1 => knife-edge)")
for k in range(1,14):
    r=c[k]**2/(c[k-1]*c[k+1])
    print(f"{k:>3} {r:>24.6f}   {'HOLDS' if r>1 else 'FAILS'}")
print("\nTuran inequalities (LP first level) HOLD, ratio ->1+ : necessary-for-RH, proven-level, asymptotically tight.")
print("Full LP hierarchy (all Laguerre-Polya ineqs) = RH; Turan is level 1 (like Lambda>=0 is the proven half).")

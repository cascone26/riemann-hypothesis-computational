#!/usr/bin/env python3
"""The single-Gaussian Weil curve D(a) = arch(a)+pole(a)−prime(a), swept over a>0.

§768/scaling finding: for a WELL-SEPARATED Gaussian-cosine basis the Weil form is ≈ diagonal, so its
positivity ⟺ every diagonal entry D(a) > 0. D(a) is the Weil functional evaluated on the single test
function φ_a(x)=e^{−x²/2}cos(ax) (Bombieri–Lagarias 'single-Gaussian' Weil positivity). D(a) > 0 ∀a is
a NECESSARY condition for RH (a single a with D(a)<0 would disprove it). Sweep a finely, find the
minimum and any near-zeros, and see whether the positivity has a floor or grazes 0.

Reuses the von-Mangoldt-correct build validated to 1e-5 vs rigorous Arb (§767 diagonal match)."""
import numpy as np, mpmath as mp
mp.mp.dps = 30
PI = np.pi; LOGPI = np.log(PI)
R = 400.0; NG = 80000                      # grid must cover peak at r=a; a up to ~360
rg = np.linspace(-R, R, NG)
print("precomputing Ω(r)...", flush=True)
OM = np.array([float(mp.re(mp.digamma(mp.mpf(1)/4 + 1j*mp.mpf(r)/2))) - LOGPI for r in rg])

def phihat(a, r):
    s2p=np.sqrt(2*PI); return s2p/2*(np.exp(-(r-a)**2/2)+np.exp(-(r+a)**2/2))
def arch(a):
    return np.trapezoid(phihat(a,rg)**2*OM, rg)/(2*PI)
def pole(a):
    s2p=np.sqrt(2*PI)
    ph=s2p/2*(np.exp(-((1j/2-a)**2)/2)+np.exp(-((1j/2+a)**2)/2))
    return 2*np.real(ph*ph)
from sympy import primerange
NMAX=500000; _p=list(primerange(2,NMAX+1)); _N=[]; _L=[]
for p in _p:
    pk=p; lp=np.log(p)
    while pk<=NMAX: _N.append(pk); _L.append(lp); pk*=p
_N=np.array(_N,float); _L=np.array(_L); _LOGN=np.log(_N); _IS=1/np.sqrt(_N)
def prime(a):
    y=_LOGN; s=2*a
    g=np.exp(-y**2/4)*(np.sqrt(PI)/2)*(np.exp(-s**2/4)+np.cos(a*y))  # diagonal g (d=0,s=2a)
    return 2*np.sum(_L*_IS*g)
def D(a): return arch(a)+pole(a)-prime(a)

# coarse sweep
aa = np.linspace(0.5, 200, 800)
DD = np.array([D(a) for a in aa])
imin = np.argmin(DD)
print(f"coarse sweep a∈[0.5,200]: min D = {DD.min():.6f} at a={aa[imin]:.3f}; "
      f"#(D≤0)={np.sum(DD<=0)}; D>0 everywhere: {np.all(DD>0)}")
# report the lowest few dips
order = np.argsort(DD)[:8]
print("  lowest 8 points (a, D):")
for k in sorted(order):
    print(f"    a={aa[k]:7.3f}  D={DD[k]:+.6f}")
# fine sweep around the global min
a0=aa[imin]; af=np.linspace(max(0.3,a0-3),a0+3,600); Df=np.array([D(a) for a in af])
print(f"fine sweep near a={a0:.2f}: min D = {Df.min():.6f} at a={af[np.argmin(Df)]:.4f}")
# small-a behavior (pole term dominant)
print("small-a: ", {round(a,1): round(D(a),4) for a in [0.5,1,2,3,5,8,12]})

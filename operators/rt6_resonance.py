#!/usr/bin/env python3
"""§452 (resonance structure of the semi-local oscillation §451): per-prime push along the FIXED critical direction
v* (bottom eigvec of the complete full). push_p = v*ᵀ prime_p v* for each prime p. Look at sign/magnitude vs log p:
is it a clean resonance (∝ Λ(p)/√p · cos(c·log p) for band center c)? That would explain the §451 oscillation as
prime-band interference — a structural handle on why positivity is global."""
import numpy as np, l2_fejer as LF
from sympy import primerange, factorint
Z=np.load('zeros150.npy'); LF.ZEROS=Z
from scipy.special import digamma
LF.RGRID=np.arange(-200.0,200.0,0.01); LF.OMEGA=np.real(digamma(0.25+0.5j*LF.RGRID))-np.log(np.pi)
BAND=(14.0,44.0); c0=(BAND[0]+BAND[1])/2; cs=LF.probes(BAND,4.0)
f=LF.build(cs); full=f['full']; NS=LF.NS
va,vv=np.linalg.eigh(full); vstar=vv[:,0]        # critical direction of the COMPLETE form
def prime_p(p):
    sc=np.zeros_like(LF.LAMW)
    for i,n in enumerate(NS):
        n=int(n); fc=factorint(n)
        if len(fc)==1 and min(fc)==p: sc[i]=1.0
    return LF.build(cs,lam_scale=sc)['prime']
primes=list(primerange(2,110))
print(f"§452 per-prime push v*ᵀ prime_p v* along critical dir. band center c0={c0}.")
print(f"{'p':>4} {'log p':>6} {'push_p':>11} {'push/(Λ/√p)':>12} {'cos(c0·logp)':>12}")
data=[]
for p in primes:
    pp=float(vstar@prime_p(p)@vstar); w=np.log(p)/np.sqrt(p)
    norm=pp/w if w>0 else 0; cc=np.cos(c0*np.log(p))
    data.append((p,pp,norm,cc))
    if p<50 or p in (53,59,71,97,109):
        print(f"{p:>4} {np.log(p):6.3f} {pp:+11.4e} {norm:+12.4f} {cc:+12.4f}")
# correlation: is normalized push ∝ cos(c0 log p)?
norms=np.array([d[2] for d in data]); ccs=np.array([d[3] for d in data])
print(f"\ncorr(normalized push, cos(c0·log p)) = {np.corrcoef(norms,ccs)[0,1]:+.3f}")
# try best-fit frequency: push_p/w ~ A cos(ω log p + φ); scan ω near band
best=(0,0)
for om in np.linspace(10,48,400):
    cc=np.cos(om*np.log(np.array(primes)))
    r=abs(np.corrcoef(norms,cc)[0,1])
    if r>best[1]: best=(om,r)
print(f"best-fit resonance frequency ω={best[0]:.2f} (corr {best[1]:.3f}); band=[{BAND[0]},{BAND[1]}], center {c0}")
print("""
reading: high corr with cos(ω log p) at ω≈band => the semi-local oscillation IS prime-band resonance: each prime
pushes ±Λ(p)/√p·cos(ω log p), and positivity requires the SIGNED SUM over all primes to be controlled — a
structured (not random) cancellation governed by the band/archimedean frequency. Low corr => the per-prime push is
not a simple single-frequency resonance (richer structure).
""")

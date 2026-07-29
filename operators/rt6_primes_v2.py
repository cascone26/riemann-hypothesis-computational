#!/usr/bin/env python3
"""§438b: recover Λ(n) from zeros with ADAPTIVE σ (scaled to resolve the local log-spacing ~1/n). σ_u = fraction
of the gap to neighboring integers. Calibrate on n=2, predict small prime-powers & composites."""
import numpy as np
from scipy.special import digamma
from scipy import integrate
from sympy import factorint
Z=np.load('zeros150.npy')
def Om(r): return np.real(digamma(0.25+0.5j*r))-np.log(np.pi)
def recover(n0):
    L=np.log(n0)
    gap=min(L-np.log(n0-1), np.log(n0+1)-L)   # local log-spacing
    sig=0.40*gap
    h=lambda r: np.exp(-sig**2*r**2/2)*np.cos(r*L)
    lhs=2*np.sum(h(Z))
    arch,_=integrate.quad(lambda r: h(r)*Om(r), -600,600, limit=900); arch/=(2*np.pi)
    pole=2*np.exp(sig**2/8.0)*np.cosh(L/2.0)
    P=(arch+pole-lhs)/2.0
    gpeak=1.0/(2*np.sqrt(2*np.pi)*sig)
    return P/gpeak*np.sqrt(n0), sig
raw2,_=recover(2); cal=np.log(2)/raw2
print(f"§438b adaptive-σ prime recovery from 150 zeros. calib const={cal:.4f}")
print(f"{'n':>3} {'σ_u':>6} {'Λ_rec':>9} {'Λ_true':>8} {'err%':>7}")
for n in [2,3,4,5,7,9,11,13,6,10,12,15]:
    rec,sig=recover(n); rec*=cal
    f=factorint(n); true=np.log(list(f.keys())[0]) if len(f)==1 else 0.0
    err=abs(rec-true)/(true+1e-9)*100 if true>0 else abs(rec)
    print(f"{n:>3} {sig:6.3f} {rec:9.4f} {true:8.4f} {err:7.1f}{'  <-composite Λ=0' if true==0 else ''}")
print("\nreading: prime-powers→log p, composites→~0 with adaptive σ => clean primes-from-zeros. else resolution-limited.")

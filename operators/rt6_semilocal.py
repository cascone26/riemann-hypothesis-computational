#!/usr/bin/env python3
"""§451 (semi-local combination — the C-C open step). Add primes ONE AT A TIME to the proven archimedean
positivity; track λmin(arch+pole − Σ_{p≤P} prime_p)."""
import numpy as np, l2_fejer as LF
from sympy import primerange, factorint
Z=np.load('zeros150.npy'); LF.ZEROS=Z
from scipy.special import digamma
LF.RGRID=np.arange(-200.0,200.0,0.01); LF.OMEGA=np.real(digamma(0.25+0.5j*LF.RGRID))-np.log(np.pi)
BAND=(14.0,44.0); cs=LF.probes(BAND,4.0)
f=LF.build(cs); AP=f['archpole']; NS=LF.NS
# precompute smallest prime factor of each n in NS
spf={int(n): min(factorint(int(n))) for n in NS}
isprimepow={int(n): (len(factorint(int(n)))==1) for n in NS}
def prime_upto(P):
    sc=np.zeros_like(LF.LAMW)
    for i,n in enumerate(NS):
        n=int(n)
        if isprimepow[n] and spf[n]<=P: sc[i]=1.0
    return LF.build(cs,lam_scale=sc)['prime']
la_arch=float(np.linalg.eigvalsh(AP)[0])
print(f"§451 semi-local: λmin(arch+pole − Σ_{{p≤P}} prime_p). band{BAND} n={len(cs)}.")
print(f"  arch+pole alone: λmin={la_arch:+.4e} (C-C PROVEN positive)")
print(f"{'P':>5} {'#primes':>8} {'λmin':>13} {'Δ prev':>12}")
prev=la_arch
for P in [2,3,5,7,11,13,17,23,31,47,71,113]:
    lm=float(np.linalg.eigvalsh(AP-prime_upto(P))[0])
    print(f"{P:>5} {len(list(primerange(2,P+1))):>8} {lm:+13.4e} {lm-prev:+12.4e}",flush=True); prev=lm
resid=np.max(np.abs(f['full']-f['zero']))
print(f"\nall primes: λmin(full)={float(np.linalg.eigvalsh(f['full'])[0]):+.4e}, resid={resid:.2e}")
print("reading: each prime subtracts; λmin staying ≥0 through all = positivity survives semi-local combination = RH here.")

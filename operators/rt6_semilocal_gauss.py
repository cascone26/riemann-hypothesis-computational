#!/usr/bin/env python3
"""§453: cross-check §451 (partial prime sums indefinite, only complete sum positive) on the GAUSSIAN instrument
(l2_edge_sweep) — the two-instrument discipline on the key new result."""
import numpy as np, l2_edge_sweep as L
from sympy import primerange, factorint
BAND=(14.0,44.0); cs=L.probes(BAND,4.0); w=3.0
f=L.build_forms(cs,w); AP=f['archpole']; NS=L.NS
isprimepow={int(n):(len(factorint(int(n)))==1) for n in NS}
spf={int(n):min(factorint(int(n))) for n in NS}
def prime_upto(P):
    sc=np.zeros_like(L.LAMW)
    for i,n in enumerate(NS):
        n=int(n)
        if isprimepow[n] and spf[n]<=P: sc[i]=1.0
    return L.build_forms(cs,w,lam_scale=sc)['prime']
la=float(np.linalg.eigvalsh(AP)[0])
print(f"§453 Gaussian instrument. arch+pole λmin={la:+.4e} (positive)")
print(f"{'P':>5} {'λmin':>13}")
for P in [2,3,5,7,11,13,17,23,31,47,71,113,199]:
    lm=float(np.linalg.eigvalsh(AP-prime_upto(P))[0])
    print(f"{P:>5} {lm:+13.4e}",flush=True)
print(f"all primes: λmin(full)={float(np.linalg.eigvalsh(f['full'])[0]):+.4e}")
print("=> if partial sums go NEGATIVE (indefinite) and only the complete sum is ≥0, §451's 'irreducibly global' holds on Gaussian too (instrument-independent).")

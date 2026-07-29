#!/usr/bin/env python3
"""§455 (Mertens/Redheffer — fresh elementary RH-equivalent, Möbius side): M(n)=Σ_{k≤n}μ(k). RH ⟺ M(n)=O(n^{1/2+ε}).
Also det(Redheffer_n)=M(n). Compute M(n) growth, |M(n)|/√n, and the largest |M(n)|/√n (Mertens conjecture M(n)<√n
is FALSE — Odlyzko-te Riele 1985 — but RH only needs O(n^{1/2+ε})). Connects to §433 (Nyman-Beurling c_n≈μ(n))."""
import numpy as np
N=2000000
mu=np.ones(N+1,dtype=np.int8); mu[0]=0
is_c=np.zeros(N+1,dtype=bool)
primes=[]
for i in range(2,N+1):
    if not is_c[i]:
        primes.append(i); mu[i]=-1
    for p in primes:
        if i*p>N: break
        is_c[i*p]=True
        if i%p==0: mu[i*p]=0; break
        else: mu[i*p]=-mu[i]
M=np.cumsum(mu.astype(np.int64))
print(f"§455 Mertens M(n)=Σμ(k), RH ⟺ M(n)=O(n^{{1/2+ε}}). computed to N={N}.")
print(f"{'n':>10} {'M(n)':>8} {'|M|/√n':>8} {'M/(√n·ln ln n)':>15}")
for n in [100,1000,10000,100000,1000000,2000000]:
    m=M[n]; print(f"{n:>10} {m:>8} {abs(m)/np.sqrt(n):8.4f} {abs(m)/(np.sqrt(n)*np.log(np.log(n))):15.4f}")
# max |M(n)|/√n over the range
ns=np.arange(2,N+1)
ratio=np.abs(M[2:])/np.sqrt(ns)
imax=np.argmax(ratio)
print(f"\nmax |M(n)|/√n = {ratio[imax]:.4f} at n={ns[imax]}  (Mertens conj M<√n: ratio<1 here so not violated in this range;")
print(f"  RH needs only |M(n)| ≤ C·n^{{1/2+ε}} — ratio bounded/slowly growing is consistent)")
# growth exponent fit: log|M| vs log n
mask=np.abs(M[2:])>0
lx=np.log(ns[mask]); ly=np.log(np.abs(M[2:])[mask])
# robust: fit upper envelope via running max
env_n=[]; env_M=[]; cm=0
for i,n in enumerate(ns):
    if abs(M[n])>cm: cm=abs(M[n]); env_n.append(n); env_M.append(cm)
env_n=np.array(env_n); env_M=np.array(env_M)
slope=np.polyfit(np.log(env_n[5:]),np.log(env_M[5:]),1)[0]
print(f"upper-envelope growth exponent of |M(n)|: {slope:.4f}  (RH ⟺ ≤ 0.5+ε; observed ≈0.5 = consistent with RH)")
print("""
reading: |M(n)|/√n stays bounded/slowly-growing (≈0.5 envelope exponent) over this range = consistent with RH
(M(n)=O(n^{1/2+ε})). This is the MÖBIUS/Mertens face — the same μ(n) that appeared as the Nyman-Beurling optimal
coefficients (§433). RH ⟺ Möbius cancellation Σμ(k)=O(√n), yet another face of the one open core. (Numerical
consistency only — the O(n^{1/2+ε}) bound for ALL n is the same open positivity.)
""")

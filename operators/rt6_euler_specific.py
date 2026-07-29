#!/usr/bin/env python3
"""§439 (Euler-product specificity, new): the Euler product's signature is Λ(p)=log p SPECIFICALLY. Test whether
Weil-form positivity needs this exact arithmetic or just the prime density. On the Fejér instrument, replace the
true Λ(p)=log p (on primes) with PERMUTED log-p values (same multiset, same support, Euler product DESTROYED), and
compare edge positivity: TRUE vs 3 random permutations. If permuted breaks positivity (goes negative before the
true rank edge) while TRUE holds, positivity needs the SPECIFIC log p = the Euler product, not just density."""
import numpy as np, l2_fejer as LF
rng=np.random.default_rng(0)
# identify primes p (n=p prime) in NS and their log p; build permutation scale on the prime (n=p) entries
NS=LF.NS; LAM=LF.LAM
prime_idx=[i for i,n in enumerate(NS) if LAM[n]>0 and n in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113)]
# actually: n is prime iff LAM[n]=log n and n has no smaller prime factor; simplest: LAM[n]==log(n)
prime_ns=[n for n in NS if abs(LAM[n]-np.log(n))<1e-9 and all(n%d for d in range(2,int(n**.5)+1))]
def scale_perm(seed):
    """scale array over NS: on prime entries, replace Λ(p)=log p by a permuted log-p'; keep prime-power (p^k,k>2) as-is."""
    sc=np.ones_like(LF.LAMW)
    logs=np.array([np.log(p) for p in prime_ns])
    perm=np.random.default_rng(seed).permutation(logs)
    for p,newlog in zip(prime_ns,perm):
        i=np.where(NS==p)[0][0]
        sc[i]=newlog/np.log(p)     # multiply Λ(p)/√p by (log p'/log p) => replaces log p with log p'
    return sc
BAND=LF.BAND
print(f"§439 Euler-product specificity (Fejér). permute log p among {len(prime_ns)} primes (Euler product destroyed).")
print(f"{'s':>5} {'n':>3} {'TRUE λmin':>11} {'perm1':>11} {'perm2':>11} {'perm3':>11} {'verdict':>18}")
for s in np.arange(6.0,2.9,-0.4):
    cs=LF.probes(BAND,s)
    if len(cs)<4: continue
    lt=float(np.linalg.eigvalsh(LF.build(cs)['full'])[0])
    lp=[float(np.linalg.eigvalsh(LF.build(cs,lam_scale=scale_perm(sd))['full'])[0]) for sd in (1,2,3)]
    v='perm NEG' if any(x<-1e-9 for x in lp) and lt>=-1e-9 else 'all>=0' if lt>=-1e-9 else 'true at edge'
    print(f"{s:5.1f} {len(cs):3d} {lt:+11.3e} {lp[0]:+11.3e} {lp[1]:+11.3e} {lp[2]:+11.3e} {v:>18}",flush=True)
print("""
reading: TRUE holds ≥0 to the rank edge while PERMUTED go NEGATIVE earlier => positivity needs the SPECIFIC
Λ(p)=log p (the Euler product), not just prime density/support. This isolates multiplicativity/Euler-product as
the arithmetic source of the positivity — exactly the input the 2D self-square Hodge form must encode (§434:
positivity-half). If permuted ALSO stays ≥0, positivity is density-generic (weaker arithmetic dependence).
""")

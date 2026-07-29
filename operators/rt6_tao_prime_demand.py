"""Tao: does the singular-series Gram's completeness margin depend on TOP-PRIME
scale (~offset s) or on the PRIMORIAL PRODUCT (~N)? Test: hold s fixed, grow the
prime set p_max in the model Gram, watch lam_min(R). If it PLATEAUS at p_max~2s
(larger primes contribute nothing), the prime-content demand is OFFSET-scale =
top-prime = fork I's favorable branch. If it keeps moving, primorial = fork II."""
import json, numpy as np
from sympy import primerange
def wf(res,p):
    k=len(res)
    if k==0: return 1.0
    nu=len(set(res))
    return (1-nu/p)/((1-1/p)**k)
def triples_at_sum(s,lo=1):
    return [(a,b,s-a-b) for a in range(lo,s) for b in range(a+1,s) if s-a-b>b]
def admissible(t):
    for p in (3,5,7):
        if len(set(x%p for x in t))==p: return False
    return True
def gram_model(T,primes):
    n=len(T); G=np.ones((n,n))
    for i in range(n):
        for j in range(n):
            X=list(T[i])+list(T[j]); v=1.0
            for p in primes: v*=wf([x%p for x in X],p)
            G[i,j]=v
    d=np.sqrt(np.diag(G))
    keep=d>1e-12
    Gk=G[np.ix_(keep,keep)]; dd=np.sqrt(np.diag(Gk))
    R=Gk/np.outer(dd,dd)
    return np.linalg.eigvalsh(R)[0], keep.sum()

for s in [42,60]:
    T=[t for t in triples_at_sum(s) if admissible(t)]
    print(f"\n=== s={s}, 2s={2*s}, {len(T)} admissible triples ===")
    print(f"{'p_max':>6} {'n_acc':>6} {'lam_min(R)':>11} {'delta':>9}")
    prev=None; rows=[]
    for pmax in [7,11,13,17,23,29,37,43,53,67,83]:
        primes=[p for p in primerange(3,pmax+1)]
        lam,nacc=gram_model(T,primes)
        delta=(lam-prev) if prev is not None else 0.0
        print(f"{pmax:>6} {nacc:>6} {lam:>11.5f} {delta:>+9.5f}")
        rows.append({"p_max":pmax,"lam_min_R":float(lam),"delta":float(delta)}); prev=lam
    json.dump({"s":s,"two_s":2*s,"ladder":rows},
              open(f"results/rt6_tao_prime_demand_s{s}.json","w"),indent=2)
print("\nPLATEAU at p_max~2s => prime demand is OFFSET-scale (top-prime, fork I favorable)")
print("wrote results/rt6_tao_prime_demand_s{42,60}.json")

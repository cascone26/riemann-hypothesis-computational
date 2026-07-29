"""Tao two-window guard (Dyson-ordered): is convergent-not-primorial N-STABLE?
Model Gram (singular series, N-independent) has a fixed all-primes limit lam*.
If the PHYSICAL lam_min(R) approaches lam* as the window grows (N, 4N, 16N) and
the required prime content does NOT recede toward ~N, the finding holds along the
coupled path. If physical falls away from lam* as N grows -> primorial worry lives."""
import json, numpy as np
from math import log
from sympy import primerange
def wf(res,p):
    k=len(res)
    if k==0: return 1.0
    nu=len(set(res)); return (1-nu/p)/((1-1/p)**k)
def triples_at_sum(s,lo=1):
    return [(a,b,s-a-b) for a in range(lo,s) for b in range(a+1,s) if s-a-b>b]
def admissible(t):
    for p in (3,5,7):
        if len(set(x%p for x in t))==p: return False
    return True
def lam_R(G):
    d=np.sqrt(np.diag(G)); keep=d>1e-12
    Gk=G[np.ix_(keep,keep)]; dd=np.sqrt(np.diag(Gk))
    return float(np.linalg.eigvalsh(Gk/np.outer(dd,dd))[0])
s=42; T=[t for t in triples_at_sum(s) if admissible(t)]; n=len(T)

# --- model limit (N-independent), plus which p_max is needed to match physical ---
def model_gram(primes):
    G=np.ones((n,n))
    for i in range(n):
        for j in range(n):
            X=list(T[i])+list(T[j]); v=1.0
            for p in primes: v*=wf([x%p for x in X],p)
            G[i,j]=v
    return G
prime_sets={pm:[p for p in primerange(3,pm+1)] for pm in [43,97,199]}
lam_model={pm:lam_R(model_gram(ps)) for pm,ps in prime_sets.items()}
print("model lam_min(R):  p<=43:%.4f  p<=97:%.4f  p<=199:%.4f  (converging limit)"%(lam_model[43],lam_model[97],lam_model[199]))

# --- physical Gram across a x4 window ladder ---
def sieve_L(M):
    spf=np.zeros(M+1,dtype=np.int64)
    for i in range(2,M+1):
        if spf[i]==0:
            for j in range(i,M+1,i):
                if spf[j]==0: spf[j]=i
    L=np.zeros(M+1)
    for m in range(2,M+1):
        p=int(spf[m]); q=m
        while q%p==0: q//=p
        if q==1: L[m]=log(p)
    return L
NMAX=12_000_000; Lam=sieve_L(NMAX+70); dL=Lam-1.0
def phys(N,chunk=300_000):
    G=np.zeros((n,n)); done=0
    while done<N:
        c=min(chunk,N-done); b=np.arange(done,done+c)
        Vc=np.array([dL[b+t[0]]*dL[b+t[1]]*dL[b+t[2]] for t in T]); G+=Vc@Vc.T; done+=c
    return G/N
print(f"\n{'N':>10} {'phys lam(R)':>11} {'gap to p<=199 model':>20}")
rows=[]
for N in [750_000,3_000_000,12_000_000]:
    lp=lam_R(phys(N)); gap=lam_model[199]-lp
    rows.append({"N":N,"phys_lam_R":lp,"gap_to_model":gap})
    print(f"{N:>10} {lp:>11.4f} {gap:>20.4f}")
gaps=[r["gap_to_model"] for r in rows]
verdict="CONVERGENT/N-STABLE (gap shrinks as N x4)" if gaps[0]>gaps[-1]>=-0.02 else "gap NOT shrinking -> primorial worry"
print(f"\ngaps: {gaps[0]:.4f} -> {gaps[1]:.4f} -> {gaps[2]:.4f}   VERDICT: {verdict}")
json.dump({"object":"Tao two-window N-stability guard on convergent-not-primorial",
 "version":"rt6-tao-twowindow-2026-07-28","s":s,"n":n,"model_limit_lam_R":lam_model,
 "physical_ladder":rows,"verdict":verdict},
 open("results/rt6_tao_twowindow_guard.json","w"),indent=2)
print("wrote results/rt6_tao_twowindow_guard.json")

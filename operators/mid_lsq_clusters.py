#!/usr/bin/env python3
"""
mid_lsq_clusters.py — Cluster scans at intermediate lsq resonances.

Context: Small-lsq (~300-400) has 80% prime stability floor → cluster collapse
unreliable at N=192. Large-lsq (~700-960) has 40% stability floor → γ₁ and γ₂
show meaningful cluster collapse.

QUESTION: At intermediate lsq (491-883), is the cluster collapse meaningful?
If the prime stability floor at lsq~500-900 is intermediate (~50%), some nearby
primes will diverge and the cluster collapse will partially discriminate.

Test clusters around:
- γ₄ resonance: lsq=491 (j=30). Nearby primes: 487, 491, 499, 503.
  Do 487, 499, 503 diverge at N=192?
- γ₇ resonance: lsq=859 (j=44). Nearby primes: 853, 857, 859, 863.
  Do 853, 857, 863 diverge at N=192?
- γ₃ resonance: lsq=883 (j=27). Nearby primes: 877, 881, 883.
  Do 877, 881 diverge at N=192?

Also: γ₆ resonance at lsq=948 (composite). Since 948 is composite, its stability
can't be from endpoint vanishing. What makes it special?
Check cluster around 948 for any competing stable values.

Also: γ₅ at j=39 (lsq=round(exp(2π×39/32.935)))
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

GAMMA5=32.935061587739189
j39=39; lsq_g5_j39=round(math.exp(TWO_PI*j39/GAMMA5))
print(f"γ₅ j=39 → lsq={lsq_g5_j39}")

def make_sieve(MAX_LSQ):
    sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
    for i in range(2,int(MAX_LSQ**0.5)+1):
        if sv[i]: sv[i*i::i]=False
    return np.where(sv)[0]

primes_520=make_sieve(520)
primes_880=make_sieve(880)
primes_960=make_sieve(960)
primes_2100=make_sieve(lsq_g5_j39+5)

def compute_ev0(lsq, N, primes):
    G=G_PER_N*N; L=2*math.log(math.sqrt(lsq)); ns=np.arange(-N,N+1); D=len(ns)
    R=TWO_PI*N/L+50.0; grid=np.linspace(-R,R,G)
    h=grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=np.real(digamma(0.25+0.5j*grid))-LOGPI
    c=TWO_PI*ns[:,None]/L; d=c-grid[None,:]
    sm=np.abs(d)<1e-12
    E=np.where(sm,np.sqrt(L)+0j,(1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch=(E*(w*Om)@E.conj().T).real/TWO_PI
    sh=np.sinh(L/4)**2; ni=ns[:,None]; nj=ns[None,:]
    B=32*L*sh*(L**2-16*math.pi**2*ni*nj)/((L**2+16*math.pi**2*ni**2)*(L**2+16*math.pi**2*nj**2))
    p_in=primes[primes<=lsq]; ll=[]; ww=[]
    for p in p_in:
        lp=math.log(p); pk=p
        while pk<=lsq:
            ll.append(math.log(pk)); ww.append(lp/math.sqrt(pk)); pk*=p
    if not ll: return 0.0
    lkpp=np.array(ll); wts_p=np.array(ww); K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,10000):
        k1=min(k0+10000,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]; Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd); tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])

results={}

def run_cluster(center, window, Nvals, primes, label="", star_list=[]):
    print(f"\n{'='*55}")
    print(f"Cluster around lsq={center} {label}")
    print(f"{'='*55}")
    stable_Nmax=[]
    for lsq in range(center-window, center+window+1):
        row=[]
        for N in Nvals:
            e=compute_ev0(lsq,N,primes)
            log_e=math.log10(abs(e)) if abs(e)>0 else -99
            row.append({"N":N,"log10":log_e})
        rd={r["N"]:r["log10"] for r in row}
        Nmax=max(Nvals); Nmin=min(Nvals)
        delta=rd[Nmax]-rd[Nmin]; stable=abs(delta)<2 and rd[Nmax]<-10
        if stable: stable_Nmax.append(lsq)
        verdict="STABLE***" if stable else "div"
        tag=" ← CENTER" if lsq==center else ""
        for s in star_list:
            if lsq==s: tag+=f" ← {s}"
        vals=" ".join(f"N{N}:{rd[N]:>6.2f}" for N in Nvals)
        print(f"  lsq={lsq:>4}:  {vals}  → {verdict}{tag}", flush=True)
        results[f"c{center}_lsq{lsq}"]={"lsq":lsq,"data":row}
    print(f"  STABLE set: {stable_Nmax}")
    if stable_Nmax==[center]:
        print(f"  *** CLUSTER COLLAPSE CONFIRMED: only {center} stable! ***")
    elif center in stable_Nmax:
        others=[x for x in stable_Nmax if x!=center]
        print(f"  {center} stable but also: {others}")
    else:
        print(f"  WARNING: {center} not the unique stable; stable={stable_Nmax}")
    return stable_Nmax

# γ₄ resonance: lsq=491 (j=30)
s491=run_cluster(491, 6, [32,64,128,192], primes_520, "(γ₄ j=30)")

# γ₇ resonance: lsq=859 (j=44)
s859=run_cluster(859, 6, [32,64,128,192], primes_880, "(γ₇ j=44)")

# γ₃ resonance: lsq=883 (j=27)
s883=run_cluster(883, 6, [32,64,128,192], primes_880, "(γ₃ j=27)")

# γ₆ resonance: lsq=948 (j=41, composite)
s948=run_cluster(948, 6, [32,64,128,192], primes_960, "(γ₆ j=41, composite!)")

# γ₅ j=39 hunt
print(f"\n=== γ₅ j=39 → lsq={lsq_g5_j39} (new territory) ===")
for N in [32,64,128,192]:
    e=compute_ev0(lsq_g5_j39,N,primes_2100)
    log_e=math.log10(abs(e)) if abs(e)>0 else -99
    results.setdefault(f"gamma5_j39_lsq{lsq_g5_j39}",[]).append({"N":N,"log10":log_e})
    print(f"  N={N:>3}: log={log_e:>7.2f}", flush=True)

print("\n\n=== SYNTHESIS ===")
print("Cluster results:")
print(f"  γ₄ (491): stable set = {s491}")
print(f"  γ₇ (859): stable set = {s859}")
print(f"  γ₃ (883): stable set = {s883}")
print(f"  γ₆ (948): stable set = {s948}")

clean_collapses=[c for c,s in [(491,s491),(859,s859),(883,s883),(948,s948)] if s==[c]]
non_unique=[(c,s) for c,s in [(491,s491),(859,s859),(883,s883),(948,s948)] if c in s and len(s)>1]
missing=[(c,s) for c,s in [(491,s491),(859,s859),(883,s883),(948,s948)] if c not in s]

if clean_collapses:
    print(f"\nCLEAN collapses (uniquely stable): {clean_collapses}")
if non_unique:
    print(f"Non-unique stable (center + others): {non_unique}")
if missing:
    print(f"WARNING — center not stable: {missing}")

os.makedirs("results",exist_ok=True)
with open("results/mid_lsq_clusters.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/mid_lsq_clusters.json")

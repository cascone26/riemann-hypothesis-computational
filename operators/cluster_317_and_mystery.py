#!/usr/bin/env python3
"""
cluster_317_and_mystery.py

Goal 1: Cluster scan around lsq=317 (γ₉ j=44) at N=32,128,192.
  IF only lsq=317 is stable in [311,323] → confirms γ₉ cluster collapse.
  This would confirm γ₉ resonance (9/10 zeros).

Goal 2: Cluster scans around mystery primes 701, 719, 727, 809.
  These were stable at N=320 but no known zero resonance was found.
  If they ALSO show cluster collapse within ±6, they might be zero resonances
  for higher zeros (γ₁₁..γ₅₀). If their neighbors are ALSO stable, they're
  generic arithmetic near-PSD, not zero resonances.

  Cluster logic: if lsq=p is a zero resonance, neighbors p±1..p±6 should
  diverge strongly by N=192. If lsq=p is generic prime near-PSD, some
  neighbors might ALSO be near-PSD (especially nearby composites on same
  "plateau" of prime structure).

Goal 3: Test 4 mystery primes at N=448 to see if global collapse continues.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

def make_sieve(MAX_LSQ):
    sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
    for i in range(2,int(MAX_LSQ**0.5)+1):
        if sv[i]: sv[i*i::i]=False
    return sv, np.where(sv)[0]

_, primes_350=make_sieve(350)
_, primes_830=make_sieve(830)

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

def run_cluster(center, window, Nvals, primes, label=""):
    print(f"\n{'='*55}")
    print(f"Cluster around lsq={center} {label} | N={Nvals}")
    print(f"{'='*55}")
    stable_at_Nmax=[]
    for lsq in range(center-window, center+window+1):
        row=[]
        for N in Nvals:
            e=compute_ev0(lsq,N,primes)
            log_e=math.log10(abs(e)) if abs(e)>0 else -99
            row.append({"N":N,"log10":log_e})
        rd={r["N"]:r["log10"] for r in row}
        Nmax=max(Nvals); Nmin=min(Nvals)
        delta=rd[Nmax]-rd[Nmin]; stable=abs(delta)<2 and rd[Nmax]<-10
        if stable: stable_at_Nmax.append(lsq)
        verdict="STABLE***" if stable else "div"
        tag=" ← CENTER" if lsq==center else ""
        vals=" ".join(f"N{N}:{rd[N]:>6.2f}" for N in Nvals)
        print(f"  lsq={lsq:>4}:  {vals}  → {verdict}{tag}", flush=True)
        results[f"cluster_{center}_lsq{lsq}"]={"lsq":lsq,"center":center,"data":row}
    print(f"  STABLE set: {stable_at_Nmax}")
    return stable_at_Nmax

# ====== CLUSTER 1: γ₉ resonance at lsq=317 ======
s317=run_cluster(317, 6, [32,128,192], primes_350, "(γ₉ j=44)")

# ====== CLUSTER 2: Mystery prime 701 ======
s701=run_cluster(701, 6, [64,128,192], primes_830, "(mystery stable)")

# ====== CLUSTER 3: Mystery prime 727 ======
s727=run_cluster(727, 6, [64,128,192], primes_830, "(mystery stable)")

# ====== CLUSTER 4: Mystery prime 809 ======
s809=run_cluster(809, 6, [64,128,192], primes_830, "(mystery stable)")

# ====== Goal 3: Mystery primes at N=448 (collapse test) ======
print("\n=== Mystery primes at N=448 (next collapse checkpoint) ===")
mystery_primes=[701, 719, 727, 809]
for p in mystery_primes:
    t0=time.time()
    e=compute_ev0(p,448,primes_830)
    dt=time.time()-t0
    log_e=math.log10(abs(e)) if abs(e)>0 else -99
    verdict="near-PSD" if log_e<-10 else "DIVERGING"
    print(f"  p={p}: N=448 log={log_e:>7.2f}  {dt:.1f}s  → {verdict}", flush=True)
    results[f"mystery_p{p}_N448"]={"p":p,"N":448,"log10":log_e}

print("\n\n=== SUMMARY ===")
print(f"Cluster 317 stable set: {s317}")
if s317==[317]:
    print("  *** γ₉ CLUSTER COLLAPSE CONFIRMED: only lsq=317 stable in window ***")
    print("  → γ₉ resonance at j=44, lsq=317 is CONFIRMED (9/10 zeros!)")
elif 317 in s317:
    print(f"  γ₉=317 stable but {len(s317)-1} neighbor(s) also stable: {[x for x in s317 if x!=317]}")
else:
    print("  WARNING: lsq=317 not stable in cluster scan — needs investigation")

print(f"\nCluster 701 stable set: {s701}")
if len(s701)==1 and 701 in s701:
    print("  701 UNIQUELY STABLE → probably a Riemann zero resonance (higher zero)")
else:
    print(f"  701 NOT uniquely stable ({s701}) → generic prime arithmetic near-PSD")

print(f"\nCluster 727 stable set: {s727}")
print(f"Cluster 809 stable set: {s809}")

os.makedirs("results",exist_ok=True)
with open("results/cluster_317_and_mystery.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/cluster_317_and_mystery.json")

#!/usr/bin/env python3
"""
cluster_highN.py — High-N tests to see if non-resonant primes diverge.

Issue found: clusters around 317 (γ₉) and 331 (γ₈) contain MULTIPLE stable primes at N=192.
The cluster collapse identifying zero resonances may require higher N.

Test 1: lsq=311, 313 (stable competitors of γ₉=317 in cluster [311,323]) at N=256,320.
  If 311 and 313 diverge at N=256 while 317 stays stable → γ₉ cluster collapse at N=256.

Test 2: lsq=337 (next prime after 331 = γ₈ resonance) at N=32,128,192,256.
  If 337 is NOT stable → γ₈ cluster collapse is within the narrow window tested.
  If 337 IS stable → need higher N to discriminate γ₈.

Test 3: lsq=719 at N=256,320,448 (mystery prime, still stable at N=448 per partial output).
  This establishes an upper bound on when non-resonant primes diverge.

Test 4: lsq=701 at N=512 (the most persistent mystery prime — stable to N=448).
  Compare to lsq=787 (γ₁ resonance, stable to N=512).

Test 5: primes [305, 360] mini scan at N=32,64,128,192.
  For the γ₈/γ₉ range: which primes are stable? If most are stable, the cluster
  collapse must be at N=256+.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

def make_sieve(MAX_LSQ):
    sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
    for i in range(2,int(MAX_LSQ**0.5)+1):
        if sv[i]: sv[i*i::i]=False
    return np.where(sv)[0]

primes_380=make_sieve(380)
primes_740=make_sieve(740)

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

# TEST 1: 311, 313 at N=256,320 (competitors of γ₉=317)
print("=== TEST 1: γ₉ competitors at high N ===")
print("Cluster [311,323]: 311,313 stable at N=192; need N=256 to see collapse?")
for lsq in [311,313,317]:
    label="γ₉ resonance" if lsq==317 else "competitor"
    print(f"\n  lsq={lsq} ({label}):", flush=True)
    row=[]
    for N in [256,320]:
        e=compute_ev0(lsq,N,primes_380)
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"log10":log_e})
        print(f"    N={N:>3}: log={log_e:>7.2f}", flush=True)
    results[f"lsq{lsq}_highN"]=row

print("\n  Summary:")
stable256=[lsq for lsq in [311,313,317] if any(r["log10"]<-10 for r in results[f"lsq{lsq}_highN"] if r["N"]==256)]
stable320=[lsq for lsq in [311,313,317] if any(r["log10"]<-10 for r in results[f"lsq{lsq}_highN"] if r["N"]==320)]
print(f"  Stable at N=256: {stable256}")
print(f"  Stable at N=320: {stable320}")
if stable320==[317]:
    print("  *** γ₉ CLUSTER COLLAPSE at N=320: only 317 stable! 9th zero confirmed. ***")
elif 317 in stable320:
    non317=[x for x in stable320 if x!=317]
    print(f"  317 stable but {non317} also stable at N=320 — need higher N")
else:
    print("  317 NOT stable at N=320 — γ₉ identification uncertain")

# TEST 2: lsq=337 (next prime beyond γ₈=331 cluster window)
print("\n\n=== TEST 2: lsq=337 (next prime after γ₈=331) ===")
for N in [32,128,192,256]:
    e=compute_ev0(337,N,primes_380)
    log_e=math.log10(abs(e)) if abs(e)>0 else -99
    results.setdefault("lsq337",[]).append({"N":N,"log10":log_e})
    print(f"  N={N:>3}: log={log_e:>7.2f}", flush=True)
rd337={r["N"]:r["log10"] for r in results["lsq337"]}
stable337=rd337.get(192,-99)<-10 and rd337.get(256,-99)<-10
print(f"  337 stable at N=192,256: {stable337}")
print(f"  (For context: 331 is γ₈ resonance, 337 is next prime)")

# TEST 3: Mini prime scan [305-360] at N=32,64,128,192
# To understand the FULL prime stability landscape at γ₈/γ₉ scale
print("\n\n=== TEST 3: Mini prime scan [305,360] at N=32..192 ===")
sv380=np.ones(361,dtype=bool); sv380[0]=sv380[1]=False
for i in range(2,19):
    if sv380[i]: sv380[i*i::i]=False
primes_in_range=[p for p in range(305,361) if p<361 and sv380[p]]
print(f"Primes in [305,360]: {primes_in_range}")
stable_in_range=[]
for p in primes_in_range:
    row=[]
    for N in [32,64,128,192]:
        e=compute_ev0(p,N,primes_380)
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"log10":log_e})
    results[f"miniscan_p{p}"]=row
    rd={r["N"]:r["log10"] for r in row}
    delta=rd[192]-rd[32]; stable=abs(delta)<2 and rd[192]<-10
    pseudo=rd.get(128,-99)<-10 and not stable
    v="STABLE***" if stable else ("PSEUDO" if pseudo else "div")
    if stable: stable_in_range.append(p)
    is_res=""
    if p==317: is_res=" ← γ₉ j=44"
    if p==331: is_res=" ← γ₈ j=40"
    print(f"  p={p:>4}: N32={rd[32]:>6.2f} N64={rd[64]:>6.2f} N128={rd[128]:>6.2f} N192={rd[192]:>6.2f}  → {v}{is_res}", flush=True)
print(f"\nStable primes in [305,360]: {stable_in_range}")
n_stable=len(stable_in_range)
n_total=len(primes_in_range)
print(f"  {n_stable}/{n_total} primes stable at N=192")
if set(stable_in_range)=={317,331}:
    print("  *** ONLY zero resonances stable! Zero resonances ARE unique at N=192 in this range! ***")
elif {317,331}.issubset(set(stable_in_range)):
    others=[p for p in stable_in_range if p not in {317,331}]
    print(f"  Zero resonances 317,331 stable but also: {others}")
else:
    print(f"  One or both resonances not stable: {[p for p in {317,331} if p not in stable_in_range]}")

os.makedirs("results",exist_ok=True)
with open("results/cluster_highN.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/cluster_highN.json")

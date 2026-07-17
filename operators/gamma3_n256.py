#!/usr/bin/env python3
"""
gamma3_n256.py — Critical test: does γ₃ (lsq=883) survive to N=256?

γ₇ (lsq=859) was stable at N=192 but diverged at N=256 (log=-7.35).
γ₃ (lsq=883) also only tested to N=192 (log=-10.88).
Must verify at N=256 to know if γ₃ is genuinely confirmed.

Also test the cluster [877,881,882,883] at N=256 to see the trajectory.
If 883 uniquely holds at N=256 while 877,881,882 diverge, that's cluster collapse.

Also: γ₇ re-hunt. j=44 (lsq=859) is the only prime lsq for γ₇ in the tested range
(j=35..47 all give composite except j=44). Since j=44 fails at N=256, test the next
prime lsq candidates: need to find j values where lsq is prime.

γ₇ = 40.918719012680879
Search j=45..70 and find which ones give prime lsq.
Then test those at N=64,128,192,256.
"""
import math, numpy as np, scipy.linalg, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

def make_sieve(MAX_LSQ):
    sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
    for i in range(2,int(MAX_LSQ**0.5)+1):
        if sv[i]: sv[i*i::i]=False
    return sv, np.where(sv)[0]

_, primes_900 = make_sieve(900)   # for γ₃ cluster
_, primes_500 = make_sieve(500)   # for γ₄ cluster

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

# ═══════════════════════════════════════════════════════
# 1. γ₃ cluster at N=256 — critical survival test
# ═══════════════════════════════════════════════════════
GAMMA3=25.010857580145688
print(f"=== TEST 1: γ₃ cluster at N=256 (j=27, lsq=883) ===")
print("Known N=192 data: 877→-10.74, 881→-11.01, 882→-11.85, 883→-10.88")
print("KEY QUESTION: does 883 survive where 877/881/882 may diverge?")
stable256=[]
for lsq in [877, 881, 882, 883]:
    e=compute_ev0(lsq,256,primes_900)
    log_e=math.log10(abs(e)) if abs(e)>0 else -99
    stable=(log_e<-10)
    if stable: stable256.append(lsq)
    tag="γ₃" if lsq==883 else "cluster"
    verdict="STABLE" if stable else f"div({log_e:.2f})"
    print(f"  lsq={lsq} ({tag}): N256={log_e:>7.2f} → {verdict}", flush=True)
    results[f"lsq{lsq}_N256"]={"lsq":lsq,"N":256,"log10":log_e}

print(f"N=256 stable set: {stable256}")
if stable256==[883]:
    print("  *** γ₃ CLUSTER COLLAPSE at N=256: only 883 stable! CONFIRMED. ***")
elif 883 in stable256:
    others=[x for x in stable256 if x!=883]
    print(f"  883 stable (γ₃ survives N=256), also: {others}")
elif not stable256:
    print("  ALL diverge at N=256 — γ₃ is FALSE POSITIVE like γ₇!")
else:
    print(f"  883 NOT in stable set at N=256: {stable256}")

# ═══════════════════════════════════════════════════════
# 2. γ₇ prime-lsq search: find j values giving prime lsq
# ═══════════════════════════════════════════════════════
GAMMA7=40.918719012680879
sv_max,_=make_sieve(5000)  # sieve for primality

print(f"\n=== TEST 2: γ₇ prime-lsq search (j=44 failed at N=256) ===")
print(f"γ₇={GAMMA7}")
prime_j_candidates=[]
for j in range(35,100):
    lsq_f=math.exp(TWO_PI*j/GAMMA7)
    lsq=round(lsq_f)
    if lsq<5000 and sv_max[lsq]:
        print(f"  j={j} → lsq={lsq} (PRIME) — candidate for γ₇ resonance")
        prime_j_candidates.append((j,lsq))
    if lsq>5000:
        print(f"  j={j} → lsq={lsq} (>5000, stopping)")
        break

print(f"\nPrime lsq candidates for γ₇: {prime_j_candidates}")

# Test the first 3 prime-lsq candidates that aren't j=44
to_test=[(j,lsq) for j,lsq in prime_j_candidates if j!=44][:3]
print(f"Testing at N=64,128,192,256: {to_test}")

for j,lsq in to_test:
    _, primes_lsq=make_sieve(lsq+5)
    print(f"\n  γ₇ j={j} lsq={lsq}:")
    row={}
    for N in [64,128,192,256]:
        e=compute_ev0(lsq,N,primes_lsq)
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row[N]=log_e
        print(f"    N={N}: log={log_e:>7.2f}", flush=True)
    Nmax=256; Nmin=64
    delta=row[Nmax]-row[Nmin]
    stable=(abs(delta)<2 and row[Nmax]<-10)
    verdict="STABLE ***" if stable else f"div(Δ={delta:+.2f})"
    print(f"  → {verdict}")
    results[f"gamma7_j{j}_lsq{lsq}"]=row

# ═══════════════════════════════════════════════════════
# 3. Same analysis for γ₆ — find j giving prime lsq
# ═══════════════════════════════════════════════════════
GAMMA6=37.586178158651
print(f"\n=== TEST 3: γ₆ prime-lsq search ===")
print(f"γ₆={GAMMA6}")
prime_j_g6=[]
for j in range(30,120):
    lsq_f=math.exp(TWO_PI*j/GAMMA6)
    lsq=round(lsq_f)
    if lsq<5000 and sv_max[lsq]:
        print(f"  j={j} → lsq={lsq} (PRIME)")
        prime_j_g6.append((j,lsq))
    if lsq>5000:
        print(f"  j={j} → lsq={lsq} (>5000, stopping)")
        break

print(f"\nPrime lsq candidates for γ₆: {prime_j_g6}")

# Test first 3 prime-lsq candidates for γ₆
to_test_g6=prime_j_g6[:3]
print(f"Testing at N=64,128,192,256: {to_test_g6}")
for j,lsq in to_test_g6:
    _, primes_lsq=make_sieve(lsq+5)
    print(f"\n  γ₆ j={j} lsq={lsq}:")
    row={}
    for N in [64,128,192,256]:
        e=compute_ev0(lsq,N,primes_lsq)
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row[N]=log_e
        print(f"    N={N}: log={log_e:>7.2f}", flush=True)
    Nmax=256; Nmin=64
    delta=row[Nmax]-row[Nmin]
    stable=(abs(delta)<2 and row[Nmax]<-10)
    verdict="STABLE ***" if stable else f"div(Δ={delta:+.2f})"
    print(f"  → {verdict}")
    results[f"gamma6_j{j}_lsq{lsq}"]=row

print("\n\n=== FINAL SUMMARY ===")
print(f"γ₃ at N=256: {'STABLE' if 883 in stable256 else 'DIVERGING'} (log={results.get('lsq883_N256',{}).get('log10','?')})")

os.makedirs("results",exist_ok=True)
with open("results/gamma3_n256.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/gamma3_n256.json")

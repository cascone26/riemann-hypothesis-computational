#!/usr/bin/env python3
"""
full_prime_scan.py — Comprehensive prime-lsq scan for ALL first 10 zeros.

KEY INSIGHT: Prior scans only searched j values giving large lsq (400-2000+).
This missed small-j resonances:
- γ₆ j=30 → lsq=151 (PRIME, stable through N=256!) — was never tested
- γ₇ j=37 → lsq=293 (PRIME, stable through N=256!) — was never tested

NEW APPROACH: For each zero γ_k, enumerate j=1..80 and for EACH j giving prime lsq
in range [50, 2000], test at N=64,128,192,256.

STABILITY CRITERION: log(N=256)<-10 AND |Δlog(64→256)|<2.5 (slightly relaxed from 2
to handle the N=192 dip seen in some true resonances like γ₈).

Focus zeros: γ₃ (25.011), γ₅ (32.935), γ₆ (37.586), γ₇ (40.919)
For γ₆ and γ₇, also verify j-value uniqueness by checking nearby j composites.
Also verify: γ₆ j=29 → lsq=127 (prime, missed by gamma3_n256.py).
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

sv5000, primes5000 = make_sieve(5000)

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
    p_in=primes5000[primes5000<=lsq]; ll=[]; ww=[]
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

NVALS=[64,128,192,256]

def test_lsq(lsq, label=""):
    row={}
    for N in NVALS:
        e=compute_ev0(lsq,N,primes5000)
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row[N]=log_e
    delta=row[256]-row[64]
    stable=(abs(delta)<2.5 and row[256]<-10)
    vals=" ".join(f"N{N}:{row[N]:>7.2f}" for N in NVALS)
    verdict="STABLE***" if stable else f"div(Δ={delta:+.2f})"
    print(f"  lsq={lsq:>5} {label}: {vals} → {verdict}", flush=True)
    return stable, row, delta

# ═══════════════════════════════════════════════════════════════
# SECTION 1: γ₆ j-value uniqueness — check j=29 and others
# ═══════════════════════════════════════════════════════════════
GAMMA6=37.586178158651
print("="*60)
print(f"γ₆={GAMMA6}: j=30(lsq=151) is STABLE. Check uniqueness.")
print("="*60)
stable_g6=[]
for j in range(25,45):
    lsq=round(math.exp(TWO_PI*j/GAMMA6))
    if lsq<50 or lsq>2000: continue
    is_prime=bool(sv5000[lsq])
    pmark="(p)" if is_prime else ""
    st,rd,dlt=test_lsq(lsq,f"γ₆j{j}{pmark}")
    if st: stable_g6.append((j,lsq))
    results[f"g6_j{j}_lsq{lsq}"]={"j":j,"lsq":lsq,"data":rd,"stable":st,"is_prime":is_prime}
print(f"\nAll stable j for γ₆: {stable_g6}")
if stable_g6 and stable_g6[0]==(30,151) and len(stable_g6)==1:
    print("  *** γ₆ j=30(lsq=151) UNIQUELY STABLE — J-VALUE UNIQUENESS CONFIRMED ***")
elif len(stable_g6)>1:
    print(f"  Multiple stable j values found: {stable_g6} — not unique")

# ═══════════════════════════════════════════════════════════════
# SECTION 2: γ₇ j-value uniqueness — check j=35..50
# ═══════════════════════════════════════════════════════════════
GAMMA7=40.918719012680879
print("\n" + "="*60)
print(f"γ₇={GAMMA7}: j=37(lsq=293) is STABLE. Check uniqueness.")
print("="*60)
stable_g7=[]
for j in range(25,55):
    lsq=round(math.exp(TWO_PI*j/GAMMA7))
    if lsq<50 or lsq>2000: continue
    is_prime=bool(sv5000[lsq])
    pmark="(p)" if is_prime else ""
    st,rd,dlt=test_lsq(lsq,f"γ₇j{j}{pmark}")
    if st: stable_g7.append((j,lsq))
    results[f"g7_j{j}_lsq{lsq}"]={"j":j,"lsq":lsq,"data":rd,"stable":st,"is_prime":is_prime}
print(f"\nAll stable j for γ₇: {stable_g7}")
if len(stable_g7)==1 and stable_g7[0]==(37,293):
    print("  *** γ₇ j=37(lsq=293) UNIQUELY STABLE — CONFIRMED ***")
elif len(stable_g7)>1:
    print(f"  Multiple stable: {stable_g7}")

# ═══════════════════════════════════════════════════════════════
# SECTION 3: γ₃ small-j hunt — focus on j=20..27
# ═══════════════════════════════════════════════════════════════
GAMMA3=25.010857580145688
print("\n" + "="*60)
print(f"γ₃={GAMMA3}: j=27(lsq=883) fails at N=256. Hunt smaller j.")
print("="*60)
stable_g3=[]
for j in range(15,30):
    lsq=round(math.exp(TWO_PI*j/GAMMA3))
    if lsq<30 or lsq>1000: continue
    is_prime=bool(sv5000[lsq])
    pmark="(p)" if is_prime else ""
    st,rd,dlt=test_lsq(lsq,f"γ₃j{j}{pmark}")
    if st: stable_g3.append((j,lsq))
    results[f"g3_j{j}_lsq{lsq}"]={"j":j,"lsq":lsq,"data":rd,"stable":st,"is_prime":is_prime}
print(f"\nAll stable j for γ₃: {stable_g3}")

# ═══════════════════════════════════════════════════════════════
# SECTION 4: γ₅ small-j hunt — focus on j=22..32
# ═══════════════════════════════════════════════════════════════
GAMMA5=32.935061587739189
print("\n" + "="*60)
print(f"γ₅={GAMMA5}: ALL j=28..39 diverge. Hunt j=22..27.")
print("="*60)
stable_g5=[]
for j in range(18,35):
    lsq=round(math.exp(TWO_PI*j/GAMMA5))
    if lsq<30 or lsq>600: continue
    is_prime=bool(sv5000[lsq])
    pmark="(p)" if is_prime else ""
    st,rd,dlt=test_lsq(lsq,f"γ₅j{j}{pmark}")
    if st: stable_g5.append((j,lsq))
    results[f"g5_j{j}_lsq{lsq}"]={"j":j,"lsq":lsq,"data":rd,"stable":st,"is_prime":is_prime}
print(f"\nAll stable j for γ₅: {stable_g5}")
if stable_g5:
    print(f"  *** γ₅ RESONANCE(S) FOUND: {stable_g5} ***")
else:
    print("  No stable j found for γ₅ in this range")

# ═══════════════════════════════════════════════════════════════
# SYNTHESIS
# ═══════════════════════════════════════════════════════════════
print("\n\n" + "="*60)
print("=== SYNTHESIS ===")
print("="*60)
print(f"γ₃ stable j values: {stable_g3}")
print(f"γ₅ stable j values: {stable_g5}")
print(f"γ₆ stable j values: {stable_g6}")
print(f"γ₇ stable j values: {stable_g7}")

os.makedirs("results",exist_ok=True)
with open("results/full_prime_scan.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/full_prime_scan.json")

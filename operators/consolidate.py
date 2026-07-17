#!/usr/bin/env python3
"""
consolidate.py — Fix issues found in mid_lsq_clusters and test at higher N.

Issues:
1. γ₃ cluster (883): sieve bug — make_sieve(880) misses prime 881 for lsq≥882.
   lsq=882 and 883 got wrong P-matrix. Rerun with make_sieve(900).

2. γ₆ at lsq=948: prior confirmation was N=32→128 only. At N=192, log=-9.44 > -10.
   Need: (a) test at N=256,320; (b) try j=40→lsq=802 and j=42→lsq=1120.

3. γ₄ cluster: 487 also stable at N=192 alongside 491. Test at N=256,320 to see if
   487 diverges (confirming γ₄=491 is the unique true resonance).

4. γ₇ cluster: 857,858,863,864 also stable at N=192 alongside 859. Test at N=256,320.
"""
import math, numpy as np, scipy.linalg, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

def make_sieve(MAX_LSQ):
    sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
    for i in range(2,int(MAX_LSQ**0.5)+1):
        if sv[i]: sv[i*i::i]=False
    return np.where(sv)[0]

primes_500  = make_sieve(500)
primes_900  = make_sieve(900)   # fix for γ₃ cluster (covers 877..889 correctly)
primes_960  = make_sieve(960)   # for γ₆ at 948
primes_880  = make_sieve(880)   # for γ₇ at 859 (covers 853..865 correctly)
primes_1200 = make_sieve(1200)  # for γ₆ j=40 (lsq=802) and j=42 (lsq=1120)

GAMMA6 = 37.586178158651
lsq_g6_j40 = round(math.exp(TWO_PI*40/GAMMA6))
lsq_g6_j42 = round(math.exp(TWO_PI*42/GAMMA6))
print(f"γ₆ j=40 → lsq={lsq_g6_j40}")
print(f"γ₆ j=42 → lsq={lsq_g6_j42}")

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

def run_lsq(lsq, Nvals, primes, label=""):
    row={}
    for N in Nvals:
        e=compute_ev0(lsq,N,primes)
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row[N]=log_e
    Nmin,Nmax=min(Nvals),max(Nvals)
    delta=row[Nmax]-row[Nmin]
    stable=abs(delta)<2 and row[Nmax]<-10
    vals=" ".join(f"N{N}:{row[N]:>7.2f}" for N in Nvals)
    verdict="STABLE" if stable else f"div(Δ={delta:+.2f},N{Nmax}={row[Nmax]:.2f})"
    print(f"  lsq={lsq:>5} {label}: {vals} → {verdict}", flush=True)
    results[f"lsq{lsq}_{label}"]={"lsq":lsq,"data":row,"stable":stable}
    return stable, row

# ═══════════════════════════════════════════════════════
# 1. FIX: γ₃ cluster with correct sieve (make_sieve(900))
# ═══════════════════════════════════════════════════════
print("\n=== FIX 1: γ₃ cluster at 883 — CORRECTED sieve (primes_900) ===")
print("Bug: prior run used primes_880, missing prime 881 for lsq≥882")
stable_883=[]
for lsq in range(877, 890):
    st,rd=run_lsq(lsq,[32,64,128,192],primes_900,
                  "γ₃center" if lsq==883 else "")
    if st: stable_883.append(lsq)
print(f"CORRECTED stable set in [877,889]: {stable_883}")
if stable_883==[883]:
    print("  *** γ₃ CLUSTER COLLAPSE CONFIRMED (correct sieve) ***")
elif 883 in stable_883:
    print(f"  883 stable but also: {[x for x in stable_883 if x!=883]}")
else:
    print(f"  883 NOT stable — γ₃ identification uncertain")

# ═══════════════════════════════════════════════════════
# 2. γ₄ cluster at N=256,320 — does 487 diverge?
# ═══════════════════════════════════════════════════════
print("\n=== TEST 2: γ₄ cluster at N=256,320 — 487 vs 491 ===")
print("N=192 stable set was [487, 491]. Need collapse to just [491].")
for lsq in [487, 491]:
    lab="competitor" if lsq==487 else "γ₄=j30"
    run_lsq(lsq,[192,256,320],primes_500,lab)

# ═══════════════════════════════════════════════════════
# 3. γ₇ cluster at N=256 — do competitors diverge?
# ═══════════════════════════════════════════════════════
print("\n=== TEST 3: γ₇ cluster at N=256 — 857/858/863/864 vs 859 ===")
print("N=192 stable set was [857,858,859,863,864]. Need collapse to [859].")
for lsq in [857, 858, 859, 863, 864]:
    lab="γ₇" if lsq==859 else "competitor"
    run_lsq(lsq,[192,256],primes_880,lab)

# ═══════════════════════════════════════════════════════
# 4. γ₆ at N=256,320 — does it recover?
# ═══════════════════════════════════════════════════════
print("\n=== TEST 4: γ₆ at N=256,320 — does lsq=948 stabilize? ===")
print("N=192 gave log=-9.44 (FAILS criterion). Need log<-10 to recover.")
run_lsq(948,[192,256,320],primes_960,"γ₆j41")

# ═══════════════════════════════════════════════════════
# 5. γ₆ alternate j values — j=40 and j=42
# ═══════════════════════════════════════════════════════
print(f"\n=== TEST 5: γ₆ alternate j — j=40→lsq={lsq_g6_j40}, j=42→lsq={lsq_g6_j42} ===")
for lsq,primes,lab in [
    (lsq_g6_j40, primes_960 if lsq_g6_j40<=960 else primes_1200, "γ₆j40"),
    (lsq_g6_j42, primes_1200, "γ₆j42"),
]:
    run_lsq(lsq,[32,64,128,192],primes,lab)

# ═══════════════════════════════════════════════════════
# 6. Also test lsq=487 as potential higher-zero resonance
# ═══════════════════════════════════════════════════════
print("\n=== TEST 6: lsq=487 — higher-zero resonance check ===")
print("Checking if 487 corresponds to any zero γ₁₁..γ₂₀ at some j.")
# For a zero γ to give lsq=487: j = γ × ln(487)/(2π) ≈ γ × 0.9842
# γ₁₁=52.97 → j=52.14 → not integer
# γ₁₂=56.45 → j=55.56 → not integer
# γ₁₃=59.35 → j=58.41 → not integer
# etc. — 487 likely arithmetic prime near-PSD
# Verify: test at N=320,448 to see if it diverges
run_lsq(487,[256,320],primes_500,"lsq487_highN")

print("\n\n=== SYNTHESIS ===")
print("γ₃ (883):", "STABLE" if 883 in stable_883 else "DIVERGING (bug confirmed or genuine)")
for key in results:
    r=results[key]; st=r.get("stable",False)
    if st:
        rd=r["data"]; Ns=sorted(rd.keys())
        print(f"  STABLE: {key} — {dict(rd)}")

os.makedirs("results",exist_ok=True)
with open("results/consolidate.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/consolidate.json")

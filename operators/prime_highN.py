#!/usr/bin/env python3
"""
prime_highN.py — Test the 11 stable primes from prime_scan.py at N=256, 320.

prime_scan showed 11/15 primes in [700,810] stable at N=192:
  701, 709, 719, 727, 739, 743, 757, 769, 773, 787, 809

Known zero resonances in this set:
  709: γ₁₀ (j=52) — CONFIRMED stable
  787: γ₁  (j=15) — CONFIRMED stable

THE CRITICAL QUESTION:
  Do all 11 primes REMAIN stable at N=256?
  Or do the 9 non-resonant primes diverge, leaving only 709 and 787?

If only 709 and 787 survive → cluster collapse operates globally across all primes.
If most survive → endpoint vanishing creates systematic prime stability;
  our operator does NOT discriminate zero resonances from arbitrary primes at N≤320.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900
MAX_LSQ=820
sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sv[i]: sv[i*i::i]=False
all_primes=np.where(sv)[0]

# 11 primes that were STABLE at N=192 in prime_scan.py
STABLE_PRIMES_N192=[701, 709, 719, 727, 739, 743, 757, 769, 773, 787, 809]
N_VALUES=[256, 320]

RESONANCES={709:"γ₁₀ j=52", 787:"γ₁ j=15"}

print("Testing stable primes from prime_scan at higher N...")
print(f"Primes: {STABLE_PRIMES_N192}")
print()

def compute_ev0(lsq, N):
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
    p_in=all_primes[all_primes<=lsq]; ll=[]; ww=[]
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
for lsq in STABLE_PRIMES_N192:
    row=[]
    for N in N_VALUES:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"log10_abs_ev0":log_e,"time_s":dt})
    results[lsq]=row

    rd={r["N"]:r["log10_abs_ev0"] for r in row}
    # Use N=192 value from prime_scan as prior (we didn't recompute it here)
    stable256=rd.get(256,-99)<-10
    stable320=rd.get(320,-99)<-10
    note=f"  ← {RESONANCES[lsq]}" if lsq in RESONANCES else ""
    status="STABLE" if stable256 and stable320 else ("STABLE_256_ONLY" if stable256 else "DIVERGING")
    print(f"p={lsq:>4}:  N256={rd.get(256,-99):>6.2f}  N320={rd.get(320,-99):>6.2f}  → {status}{note}", flush=True)

print()
print("=== COLLAPSE ANALYSIS ===")
surviving_256=[p for p in STABLE_PRIMES_N192 if results[p][0]["log10_abs_ev0"]<-10]
surviving_320=[p for p in STABLE_PRIMES_N192 if results[p][1]["log10_abs_ev0"]<-10]
print(f"Stable at N=192: {len(STABLE_PRIMES_N192)} primes: {STABLE_PRIMES_N192}")
print(f"Stable at N=256: {len(surviving_256)} primes: {surviving_256}")
print(f"Stable at N=320: {len(surviving_320)} primes: {surviving_320}")
print()
print("Key verdict:")
if surviving_320==[709,787] or set(surviving_320)=={709,787}:
    print("*** PERFECT COLLAPSE: only zero resonances (709,787) survive to N=320 ***")
    print("Endpoint vanishing + zero resonance = unique stable fixed point")
elif len(surviving_320)<=4:
    print(f"PARTIAL COLLAPSE: {len(surviving_320)} primes survive to N=320 (vs 11 at N=192)")
    print(f"Surviving: {surviving_320}")
    non_res=[p for p in surviving_320 if p not in RESONANCES]
    if non_res:
        print(f"Non-resonant survivors: {non_res} — need further investigation")
else:
    print(f"WEAK COLLAPSE: {len(surviving_320)}/{len(STABLE_PRIMES_N192)} primes still stable at N=320")
    print("Prime near-PSD persists — zero resonance NOT uniquely selected at N=320")

os.makedirs("results",exist_ok=True)
with open("results/prime_highN.json","w") as f:
    json.dump({str(k):v for k,v in results.items()},f,indent=2)
print("\nSaved → results/prime_highN.json")

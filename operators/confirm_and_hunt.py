#!/usr/bin/env python3
"""
confirm_and_hunt.py — Three goals:

1. CONFIRM γ₈ j=40 (lsq=331) stability at N=256,320.
   small_j_hunt showed: N32=-12.0 N64=-11.3 N96=-10.9 N128=-10.7 N192=-10.7 → STABLE
   Need high-N to distinguish true stable from pseudo.

2. CHECK γ₉ j=44 (lsq=317) at N=256.
   small_j_hunt showed: N32=-13.1 N64=-12.4 N128=-11.9 N192=-10.7 → PSEUDO (Δlog=2.4)
   If N=256 gives log<-10 AND Δlog from N=64..256 < 2, it's genuinely stable.

3. HUNT γ₅ at high lsq (j=37..40 → lsq=1169..1715).
   All j<37 for γ₅ tested — all diverge.
   j=37: lsq=1169, j=38: lsq=1417 (needs MAX_LSQ=1500), j=39: lsq=1715 (skip, too large)
   Also: j=36 tested (lsq=961) — diverging. j=37 is the new frontier.

Also: cluster check of lsq around 331 to see if γ₈ collapse occurs.
   Test lsq=[326..336] at N=64,128,192 (small cluster scan for γ₈).
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

all_primes_400=make_sieve(400)
all_primes_1200=make_sieve(1200)
all_primes_1500=make_sieve(1500)

ZEROS={5:32.935061587739189, 8:43.327073280914999, 9:48.005150881167159}

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

# ====== PART 1: CONFIRM γ₈ j=40 lsq=331 ======
print("=== PART 1: Confirm γ₈ j=40 lsq=331 at N=256,320 ===")
print("Prior: N32=-12.0 N64=-11.3 N96=-10.9 N128=-10.7 N192=-10.7")
for N in [256,320]:
    t0=time.time()
    e=compute_ev0(331,N,all_primes_400)
    dt=time.time()-t0
    log_e=math.log10(abs(e)) if abs(e)>0 else -99
    print(f"  N={N:>3}: log={log_e:>7.2f}  {dt:.1f}s", flush=True)
    results.setdefault("gamma8_j40_lsq331",[]).append({"N":N,"log10":log_e})

# ====== PART 2: γ₉ j=44 lsq=317 at N=256 ======
print("\n=== PART 2: γ₉ j=44 lsq=317 at N=256 ===")
print("Prior: N32=-13.1 N64=-12.4 N128=-11.9 N192=-10.7  Δlog=+2.4 (PSEUDO)")
for N in [256,320]:
    t0=time.time()
    e=compute_ev0(317,N,all_primes_400)
    dt=time.time()-t0
    log_e=math.log10(abs(e)) if abs(e)>0 else -99
    print(f"  N={N:>3}: log={log_e:>7.2f}  {dt:.1f}s", flush=True)
    results.setdefault("gamma9_j44_lsq317",[]).append({"N":N,"log10":log_e})

# ====== PART 3: γ₅ j=37 lsq=1169 ======
GAMMA5=ZEROS[5]; j37=37; lsq37=round(math.exp(TWO_PI*j37/GAMMA5))
print(f"\n=== PART 3a: γ₅ j={j37} → lsq={lsq37} [hunt] ===")
print("Prior for γ₅: all j=28..36 (lsq=209..961) DIVERGE.")
for N in [32,64,128,192]:
    t0=time.time()
    e=compute_ev0(lsq37,N,all_primes_1200)
    dt=time.time()-t0
    log_e=math.log10(abs(e)) if abs(e)>0 else -99
    print(f"  N={N:>3}: log={log_e:>7.2f}  {dt:.1f}s", flush=True)
    results.setdefault(f"gamma5_j{j37}_lsq{lsq37}",[]).append({"N":N,"log10":log_e})

j38=38; lsq38=round(math.exp(TWO_PI*j38/GAMMA5))
print(f"\n=== PART 3b: γ₅ j={j38} → lsq={lsq38} [hunt] ===")
if lsq38<=1500:
    for N in [32,64,128,192]:
        t0=time.time()
        e=compute_ev0(lsq38,N,all_primes_1500)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        print(f"  N={N:>3}: log={log_e:>7.2f}  {dt:.1f}s", flush=True)
        results.setdefault(f"gamma5_j{j38}_lsq{lsq38}",[]).append({"N":N,"log10":log_e})
else:
    print(f"  lsq={lsq38} > 1500, skipping")

# ====== PART 4: Mini cluster scan around lsq=331 (γ₈ resonance) ======
print("\n=== PART 4: Mini cluster around lsq=331 (γ₈ j=40) at N=32,128,192 ===")
for lsq in range(326,337):
    row=[]
    for N in [32,128,192]:
        t0=time.time()
        e=compute_ev0(lsq,N,all_primes_400)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"log10":log_e})
    rd={r["N"]:r["log10"] for r in row}
    delta=rd[192]-rd[32]; stable=abs(delta)<2 and rd[192]<-10
    verdict="STABLE***" if stable else "diverging"
    print(f"  lsq={lsq:>4}:  N32={rd[32]:>6.2f}  N128={rd[128]:>6.2f}  N192={rd[192]:>6.2f}  → {verdict}", flush=True)
    results[f"cluster331_lsq{lsq}"]=row

print("\n\n=== FINAL SUMMARY ===")
for key in ["gamma8_j40_lsq331","gamma9_j44_lsq317"]:
    if key not in results: continue
    rd={r["N"]:r["log10"] for r in results[key]}
    ns=sorted(rd.keys()); delta=rd[ns[-1]]-rd[ns[0]]
    maxlog=max(rd.values()); verdict="STABLE" if abs(delta)<2 and maxlog<-10 else "DIVERGING"
    print(f"  {key}: {dict(rd)} → {verdict}")

for j,lsq in [(j37,lsq37),(j38,lsq38)]:
    key=f"gamma5_j{j}_lsq{lsq}"
    if key not in results: continue
    rd={r["N"]:r["log10"] for r in results[key]}
    ns=sorted(rd.keys()); delta=rd[ns[-1]]-rd[ns[0]]
    maxlog=max(rd.values()); stable=abs(delta)<2 and maxlog<-10
    verdict="STABLE ***" if stable else "diverging"
    print(f"  γ₅ j={j} lsq={lsq}: {dict(rd)} → {verdict}")

os.makedirs("results",exist_ok=True)
with open("results/confirm_and_hunt.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/confirm_and_hunt.json")

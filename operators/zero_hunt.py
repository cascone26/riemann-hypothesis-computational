#!/usr/bin/env python3
"""
zero_hunt.py — Systematic search for stable j per zero.

Strategy: for each zero, scan all j giving lsq in [400, 1500].
Quick filter: compute N=32 and N=128 only. Flag where |Δlog|<2.5.
Then confirm flagged candidates with full N=32..192.

Targets: γ₅, γ₈, γ₉, γ₁₀ (stable j not yet found).
Also: γ₄ j=30 confirmation extended to N=256.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

ZEROS_HUNT=[
    (4,  30.424876125859513),  # confirm j=30/lsq=491; already found
    (5,  32.935061587739189),  # no stable j found yet
    (8,  43.327073280914999),  # no stable j found yet
    (9,  48.005150881167159),  # no stable j found yet
    (10, 49.773832477672302),  # no stable j found yet
]

MAX_LSQ=1600
sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sv[i]: sv[i*i::i]=False
primes=np.where(sv)[0]

def compute_ev0(lsq, N, G=None):
    if G is None: G=G_PER_N*N
    L=2*math.log(math.sqrt(lsq)); ns=np.arange(-N,N+1); D=len(ns)
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
    S=np.zeros(D); Sd=np.zeros(D); K_BATCH=10000
    for k0 in range(0,K,K_BATCH):
        k1=min(k0+K_BATCH,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]; Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd); tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])

LSQ_MIN=400; LSQ_MAX=1500
results={}
candidates={}  # flagged as potentially stable

for k,gamma in ZEROS_HUNT:
    print(f"\n{'='*60}")
    print(f"HUNTING γ_{k}={gamma:.6f} — scanning j for lsq in [{LSQ_MIN},{LSQ_MAX}]")
    print(f"{'='*60}")

    # Find all j in range
    js_in_range=[]
    for j in range(1,200):
        raw=math.exp(TWO_PI*j/gamma)
        lsq=round(raw)
        if lsq<LSQ_MIN: continue
        if lsq>LSQ_MAX: break
        js_in_range.append((j,lsq))

    print(f"Candidates: {[(j,lsq) for j,lsq in js_in_range]}")
    print()

    zero_results=[]
    for j,lsq in js_in_range:
        t0=time.time()
        e32=compute_ev0(lsq,32)
        e128=compute_ev0(lsq,128)
        dt=time.time()-t0
        log32=math.log10(abs(e32)) if abs(e32)>0 else -99
        log128=math.log10(abs(e128)) if abs(e128)>0 else -99
        delta=log128-log32
        stable=abs(delta)<2.5 and log128<-10
        marker=" *** CANDIDATE ***" if stable else ""
        print(f"  j={j:>3}  lsq={lsq:>5}  log32={log32:>7.2f}  log128={log128:>7.2f}  Δ={delta:>+6.2f}  {dt:.1f}s{marker}", flush=True)
        row={"j":j,"lsq":lsq,"log32":log32,"log128":log128,"delta":delta,"stable_flag":stable}
        zero_results.append(row)
        if stable:
            candidates.setdefault(k,[]).append((j,lsq,delta))

    results[k]={"gamma":gamma,"scan":zero_results}

print("\n\n=== CANDIDATE SUMMARY ===")
if not candidates:
    print("No stable candidates found in any zero.")
else:
    for k,cands in candidates.items():
        gamma=dict(ZEROS_HUNT)[k]
        print(f"\nγ_{k}={gamma:.6f}:")
        for j,lsq,delta in cands:
            print(f"  j={j}  lsq={lsq}  Δlog={delta:+.2f}  → STABLE CANDIDATE")

# Full N-scaling for confirmed candidates
if candidates:
    print("\n\n=== FULL N-SCALING FOR CANDIDATES ===")
    FULL_N=[32,64,96,128,192]
    full_results={}
    for k,cands in candidates.items():
        gamma=dict(ZEROS_HUNT)[k]
        for j,lsq,_ in cands:
            key=f"g{k}_j{j}_lsq{lsq}"
            print(f"\n--- γ_{k}={gamma:.3f}  j={j}  lsq={lsq} ---")
            row=[]
            for N in FULL_N:
                t0=time.time()
                e=compute_ev0(lsq,N)
                dt=time.time()-t0
                log_e=math.log10(abs(e)) if abs(e)>0 else -99
                row.append({"N":N,"log10_abs_ev0":log_e,"time_s":dt})
                print(f"  N={N:>3}: log={log_e:>7.2f}  {dt:.1f}s", flush=True)
            full_results[key]={"k":k,"gamma":gamma,"j":j,"lsq":lsq,"full_N":row}

    results["full_N_candidates"]=full_results

    print("\n=== FINAL VERDICT ===")
    for key,v in full_results.items():
        rd={r["N"]:r["log10_abs_ev0"] for r in v["full_N"]}
        v32=rd.get(32,99); v192=rd.get(192,99)
        delta=v192-v32
        verdict="STABLE" if abs(delta)<2 else "DIVERGING"
        print(f"  γ_{v['k']}  j={v['j']}  lsq={v['lsq']}:  Δlog={delta:+.2f}  → {verdict}")

os.makedirs("results",exist_ok=True)
with open("results/zero_hunt.json","w") as f:
    json.dump({str(k):v for k,v in results.items()},f,indent=2)
print("\nSaved → results/zero_hunt.json")

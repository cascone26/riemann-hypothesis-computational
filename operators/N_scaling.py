#!/usr/bin/env python3
"""
N_scaling.py — ev0(τ_N(lsq)) as a function of N.

Key question for RH: as N → ∞, does ev0(τ_N(lsq)) → 0?

If yes: the infinite-dimensional operator is PSD → Weil criterion holds → RH
If no: ev0 → some negative constant → Weil criterion fails at this scale → ¬RH

Tests at N = 8, 16, 32, 64, 96, 128 for selected lsq values:
  - lsq=787 (γ₁ resonance, j=15)
  - lsq=967 (γ₂ resonance, j=23)
  - lsq=926 (j=15 for fake γ=13.8, known "non-resonant")
  - lsq=500 (control: non-cluster value)
"""
import math, numpy as np, scipy.linalg, time
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900  # Grid points per N unit (constant resolution ratio)

TEST_LSQS=[
    (787, "γ₁ resonance (j=15)"),
    (967, "γ₂ resonance (j=23)"),
    (926, "fake resonance (j=15, γ=13.8)"),
    (779, "between-cluster (lsq=779, non-near-PSD)"),
    (500, "control (lsq=500)"),
]
N_VALUES=[8, 16, 32, 64, 96, 128]

MAX_LSQ=1000
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
    ll=[]; ww=[]
    for p in primes:
        if p>lsq: break
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

print("N-SCALING STUDY: ev0(τ_N(lsq)) vs N")
print(f"Grid points = {G_PER_N} × N (constant resolution ratio)")
print(f"N values: {N_VALUES}")
print()

results={}
for lsq,label in TEST_LSQS:
    print(f"\n=== lsq={lsq} ({label}) ===")
    print(f"{'N':>5}  {'G':>7}  {'ev0':>14}  {'log|ev0|':>10}  {'time':>7}")
    print("-"*50)
    row=[]
    for N in N_VALUES:
        G=G_PER_N*N
        t0=time.time()
        e=compute_ev0(lsq,N,G)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if e!=0 else -99
        row.append({"N":N,"G":G,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
        print(f"{N:>5}  {G:>7}  {e:>14.6e}  {log_e:>10.2f}  {dt:>6.1f}s", flush=True)
    results[lsq]={"label":label,"data":row}

print("\n\n=== SUMMARY: ev0 trend with N ===")
print(f"{'lsq':>6}  {'label':<40}", end="")
for N in N_VALUES: print(f"  N={N:>3}", end="")
print()
print("-"*100)
for lsq,label in TEST_LSQS:
    print(f"{lsq:>6}  {label:<40}", end="")
    for r in results[lsq]["data"]:
        print(f"  {r['log10_abs_ev0']:>5.2f}", end="")
    print()

print()
print("Interpretation:")
print("  If ev0 → 0 (log|ev0| → -∞) as N grows: operator is tending to PSD (consistent with RH)")
print("  If ev0 → const < 0 (log|ev0| → constant) as N grows: operator has fixed negative eigenvalue")

import json, os
os.makedirs("results",exist_ok=True)
with open("results/N_scaling.json","w") as f:
    json.dump({"G_per_N":G_PER_N,"N_values":N_VALUES,
               "test_lsqs":[{"lsq":l,"label":lb} for l,lb in TEST_LSQS],
               "results":{str(lsq):v for lsq,v in results.items()}},f,indent=2)
print(f"\nSaved → results/N_scaling.json")

#!/usr/bin/env python3
"""
Full γ scan: sweep γ from 13.0 to 22.0, step 0.01.
Adaptive j: target lsq≈787 (same prime set as our j=15 test for γ₁).

For each γ: j_k = max(1, round(γ × log(787) / (2π)))
            lsq_k = round(exp(2π×j_k/γ))  ← actual lsq used

Expected:
  γ=14.1347...: ev0 ≈ machine precision (γ₁, first Riemann zero)
  γ=21.0220...: ev0 ≈ machine precision (γ₂, second Riemann zero)
  All others:   ev0 << 0

Outputs:
  results/gamma_scan_full.out   (text table)
  results/gamma_scan_full.json  (data for plotting)
"""
import math, numpy as np, scipy.linalg, json, sys, time
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693; GAMMA2=21.022039638771555
L_TARGET=math.log(787); G=60000; N=64; K_BATCH=10000

# Sieve once for all lsq ≤ 1500 (covers full gamma range)
MAX_LSQ=1500
print(f"Sieving to {MAX_LSQ}...", flush=True)
sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sv[i]: sv[i*i::i]=False
# Precompute prime powers for each lsq value we'll need
primes_up_to=np.where(sv)[0]
print(f"  {len(primes_up_to)} primes ≤ {MAX_LSQ}", flush=True)

# Cache prime data per lsq
pd_cache={}
def get_pd(lsq):
    lsq=min(lsq,MAX_LSQ)
    if lsq not in pd_cache:
        ll,ww=[],[]
        for p in primes_up_to:
            if p>lsq: break
            lp=math.log(p); pk=p
            while pk<=lsq:
                ll.append(math.log(pk)); ww.append(lp/math.sqrt(pk)); pk*=p
        pd_cache[lsq]=(np.array(ll),np.array(ww))
    return pd_cache[lsq]

def build_ev0(gamma, N=64):
    j=max(1,round(gamma*L_TARGET/TWO_PI))
    lsq=min(MAX_LSQ,round(math.exp(TWO_PI*j/gamma)))
    if lsq<2: return float('nan'),j,lsq
    L=2*math.log(math.sqrt(lsq))
    ns=np.arange(-N,N+1); D=len(ns)
    R=TWO_PI*N/L+50.0; grid=np.linspace(-R,R,G)
    h=grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=np.real(digamma(0.25+0.5j*grid))-LOGPI
    c=TWO_PI*ns[:,None]/L; d=c-grid[None,:]
    sm=np.abs(d)<1e-12
    E=np.where(sm,np.sqrt(L)+0j,(1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch=(E*(w*Om)@E.conj().T).real/TWO_PI
    sh=np.sinh(L/4)**2; ni=ns[:,None]; nj=ns[None,:]
    B=32*L*sh*(L**2-16*math.pi**2*ni*nj)/((L**2+16*math.pi**2*ni**2)*(L**2+16*math.pi**2*nj**2))
    lkpp,wts_p=get_pd(lsq); K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,K_BATCH):
        k1=min(k0+K_BATCH,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]
        Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd)
    tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0]),j,lsq

# Gamma range
gammas=np.arange(13.0,22.01,0.01)
print(f"\nScanning {len(gammas)} γ values from {gammas[0]:.2f} to {gammas[-1]:.2f}")
print(f"Adaptive j targeting lsq≈787 (L≈{L_TARGET:.3f}), N={N}, G={G}")
print(f"\n{'γ':>8}  {'j':>4}  {'lsq':>6}  {'ev0':>14}  {'log|ev0|':>10}  flag")
print("-"*60, flush=True)

results=[]
t0=time.time()
for gam in gammas:
    ev,j,lsq=build_ev0(gam)
    log_ev=math.log10(abs(ev)) if ev!=0 else -99
    flag=""
    if abs(gam-GAMMA1)<0.015: flag="<-- γ₁ (REAL ZERO)"
    elif abs(gam-GAMMA2)<0.015: flag="<-- γ₂ (REAL ZERO)"
    elif log_ev < -10: flag="!!! NEAR-PSD !!!"  # ev0 close to zero = potentially a zero
    always_print=(log_ev<-10 or abs(gam-GAMMA1)<0.1 or abs(gam-GAMMA2)<0.1)
    if flag or always_print:
        print(f"{gam:>8.4f}  {j:>4}  {lsq:>6}  {ev:>14.4e}  {log_ev:>10.2f}  {flag}")
        sys.stdout.flush()
    results.append({"gamma":round(float(gam),4),"j":j,"lsq":int(lsq),"ev0":ev,"log10_abs_ev0":log_ev})

dt=time.time()-t0
print(f"\nDone in {dt:.1f}s ({dt/len(gammas)*1000:.1f}ms/point)")

# Find near-PSD values (potential zeros found)
near_psd=[(r["gamma"],r["ev0"]) for r in results if r["log10_abs_ev0"]<-10]
print(f"\nNear-PSD values (log|ev0|>-7): {len(near_psd)}")
for g,e in near_psd:
    print(f"  γ={g:.4f}: ev0={e:.4e}")

import os; os.makedirs("results",exist_ok=True)
with open("results/gamma_scan_full.json","w") as f:
    json.dump({"N":N,"G":G,"L_target":L_TARGET,"results":results},f)
print("\nSaved → results/gamma_scan_full.json")

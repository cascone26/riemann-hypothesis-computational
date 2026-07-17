#!/usr/bin/env python3
"""
lsq_landscape.py — ev0(τ_N(lsq)) for lsq = 100 to MAX_LSQ.

This answers the key question: what is the full landscape of the minimum
eigenvalue of the CvS operator as a function of scale lsq?

Expected output: ev0 oscillates near 0 (near-PSD clusters) with periodic
gaps of more-negative values. Under RH, the GLOBAL infimum of ev0(lsq)
over all lsq should → 0 as N → ∞.

This is the fundamental observable: not "does γ give near-PSD" but
"how negative does ev0(lsq) get across all scales?"

Output: results/lsq_landscape.json, results/lsq_landscape.out
"""
import math, numpy as np, scipy.linalg, json, sys, time
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
N=32   # Use N=32 for speed (survey); final paper uses N=64
G=30000  # Half grid points for speed

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--lsq_min",type=int,default=100)
parser.add_argument("--lsq_max",type=int,default=2000)
parser.add_argument("--step",type=int,default=1)
parser.add_argument("--N",type=int,default=32)
parser.add_argument("--G",type=int,default=30000)
args=parser.parse_args()
N=args.N; G=args.G

# Sieve
print(f"Sieving to {args.lsq_max}...", flush=True)
sv=np.ones(args.lsq_max+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(args.lsq_max**0.5)+1):
    if sv[i]: sv[i*i::i]=False
primes=np.where(sv)[0]
print(f"  {len(primes)} primes")

def ev0(lsq):
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
    # Build prime data for this lsq
    ll=[]; ww=[]
    for p in primes:
        if p>lsq: break
        lp=math.log(p); pk=p
        while pk<=lsq:
            ll.append(math.log(pk)); ww.append(lp/math.sqrt(pk)); pk*=p
    if not ll: return 0.0
    lkpp=np.array(ll); wts_p=np.array(ww); K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D)
    K_BATCH=10000
    for k0 in range(0,K,K_BATCH):
        k1=min(k0+K_BATCH,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]; Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd); tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])

lsq_vals=list(range(args.lsq_min, args.lsq_max+1, args.step))
print(f"Computing ev0 for {len(lsq_vals)} lsq values ({args.lsq_min}..{args.lsq_max}, step={args.step})")
print(f"N={N}, G={G}")
print()

import os; os.makedirs("results",exist_ok=True)
outf=open("results/lsq_landscape.out","w")
outf.write(f"# lsq_landscape: N={N}, G={G}\n")
outf.write(f"# lsq  ev0  log10|ev0|  is_prime\n")

results=[]
t0=time.time()
NEAR_PSD_THRESH=-10.0

near_psd_count=0
strongly_neg_count=0  # log|ev0| > -5 (ev0 < -1e-5)
total=len(lsq_vals)

print(f"{'lsq':>6}  {'ev0':>14}  {'log|ev0|':>10}  prime  flag")
print("-"*55, flush=True)

for i,lsq in enumerate(lsq_vals):
    e=ev0(lsq)
    log_e=math.log10(abs(e)) if e!=0 else -99
    ip=bool(sv[lsq])
    flag=""
    if abs(e)>0 and log_e<NEAR_PSD_THRESH:
        flag="NEAR-PSD"; near_psd_count+=1
    elif abs(e)>0 and log_e>-5:
        flag="NEGATIVE"; strongly_neg_count+=1
    results.append({"lsq":lsq,"ev0":e,"log10_abs_ev0":log_e,"is_prime":ip})
    outf.write(f"{lsq}  {e:.6e}  {log_e:.2f}  {ip}\n")
    # Print every 50th or notable
    if flag or i%50==0:
        print(f"{lsq:>6}  {e:>14.6e}  {log_e:>10.2f}  {ip!s:>5}  {flag}", flush=True)
    if (i+1)%200==0:
        el=time.time()-t0; rate=(i+1)/el
        pct=100*(i+1)/total
        print(f"  [{i+1}/{total}] {pct:.1f}% done, {(total-i-1)/rate:.0f}s remaining", flush=True)

outf.close()
dt=time.time()-t0

# Summary
log_evs=[r["log10_abs_ev0"] for r in results]
print(f"\n{'='*60}")
print(f"DONE in {dt:.1f}s")
print(f"lsq range: {args.lsq_min}..{args.lsq_max}")
print(f"Near-PSD (log|ev0|<-10): {near_psd_count}/{total} = {100*near_psd_count/total:.1f}%")
print(f"Strongly negative (log|ev0|>-5): {strongly_neg_count}/{total} = {100*strongly_neg_count/total:.1f}%")
print(f"Worst ev0: {10**max(log_evs):.4e} (log={max(log_evs):.2f}) at lsq={results[log_evs.index(max(log_evs))]['lsq']}")
print(f"Best ev0:  {10**min(log_evs):.4e} (log={min(log_evs):.2f}) at lsq={results[log_evs.index(min(log_evs))]['lsq']}")
print(f"Median log|ev0|: {np.median(log_evs):.2f}")
print(f"Mean log|ev0|:   {np.mean(log_evs):.2f}")

with open("results/lsq_landscape.json","w") as f:
    json.dump({"N":N,"G":G,"lsq_min":args.lsq_min,"lsq_max":args.lsq_max,
               "step":args.step,"near_psd_count":near_psd_count,
               "total":total,"near_psd_fraction":near_psd_count/total,
               "worst_log_ev0":max(log_evs),"best_log_ev0":min(log_evs),
               "results":results},f)
print(f"Saved → results/lsq_landscape.json")

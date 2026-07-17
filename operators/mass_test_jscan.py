#!/usr/bin/env python3
"""
Mass test v2: J-scanning per zero.

For each zero γ_k, scan j from 8 to 35 and find the best (least negative) ev0.
Report: best_j, best_lsq, best_ev0 for each zero.

Expected: ALL real zeros have some j giving ev0 near machine precision.
Expected: ALL fake midpoints have ev0 < -8 at ALL j values (no near-PSD).

This tests the strongest form of the Weil criterion:
  ∀k ∃j: ev0(τ_N(lsq=round(exp(2πj/γ_k)))) ≈ 0

Multiprocessed on all available cores.
Outputs:
  results/mass_test_jscan.json
  results/mass_test_jscan.out
"""
import math, numpy as np, scipy.linalg, json, sys, time, os
import multiprocessing as mp
from scipy.special import digamma
import argparse

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G=60000; N=64; K_BATCH=10000
J_RANGE=list(range(8,36))  # j values to scan

parser=argparse.ArgumentParser()
parser.add_argument("--nzeros",type=int,default=1000)
parser.add_argument("--ncores",type=int,default=os.cpu_count())
parser.add_argument("--jmin",type=int,default=8)
parser.add_argument("--jmax",type=int,default=35)
args=parser.parse_args()
J_RANGE=list(range(args.jmin,args.jmax+1))

# Load zeros
ZEROS_FILE=os.path.expanduser("~/projects/riemann/data/zeros_100k.txt")
with open(ZEROS_FILE) as f:
    all_zeros=[float(l.strip()) for l in f if l.strip()]
zeros=sorted(all_zeros[:args.nzeros])
fakes=[(zeros[i]+zeros[i+1])/2 for i in range(len(zeros)-1)]
print(f"Loaded {len(zeros)} zeros, {len(fakes)} fakes")

# Compute all lsq values needed
all_lsqs=set()
for g in zeros+fakes:
    for j in J_RANGE:
        lsq=round(math.exp(TWO_PI*j/g))
        if 4<=lsq<=4000:
            all_lsqs.add(lsq)
MAX_LSQ=max(all_lsqs)
print(f"Unique lsq values: {len(all_lsqs)}, MAX={MAX_LSQ}")

# Sieve
print(f"Sieving to {MAX_LSQ}...", flush=True)
t0=time.time()
sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sv[i]: sv[i*i::i]=False
primes=np.where(sv)[0]
print(f"  {len(primes)} primes in {time.time()-t0:.1f}s")

# Build prime data for all lsq values
print(f"Building prime data for {len(all_lsqs)} lsq values...", flush=True)
t1=time.time()
pd_map={}
for lsq in sorted(all_lsqs):
    ll,ww=[],[]
    for p in primes:
        if p>lsq: break
        lp=math.log(p); pk=p
        while pk<=lsq:
            ll.append(math.log(pk)); ww.append(lp/math.sqrt(pk)); pk*=p
    pd_map[lsq]=(np.array(ll),np.array(ww))
print(f"  Done in {time.time()-t1:.1f}s")

def ev0_at(lsq):
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
    lkpp,wts_p=pd_map[lsq]; K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,K_BATCH):
        k1=min(k0+K_BATCH,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]; Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd); tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])

def worker_jscan(args):
    """For one gamma, scan all j and return best ev0."""
    gamma,is_real=args
    results={}
    for j in J_RANGE:
        lsq=round(math.exp(TWO_PI*j/gamma))
        if lsq<4 or lsq>MAX_LSQ: continue
        if lsq not in pd_map: continue
        try:
            e=ev0_at(lsq)
            results[j]=(lsq,e)
        except: pass
    if not results: return (gamma,is_real,None,None,float('-inf'))
    best_j=max(results,key=lambda j_: results[j_][1])
    best_lsq,best_ev=results[best_j]
    return (gamma,is_real,best_j,best_lsq,best_ev)

# Build tasks
tasks=[(g,True) for g in zeros]+[(g,False) for g in fakes]
print(f"\nTotal tasks: {len(tasks)} ({len(zeros)} zeros + {len(fakes)} fakes)")
print(f"J range: {J_RANGE[0]}..{J_RANGE[-1]}, N={N}, G={G}")
print(f"Using {args.ncores} cores", flush=True)

os.makedirs("results",exist_ok=True)
outf=open("results/mass_test_jscan.out","w")
outf.write(f"# mass_test_jscan: j={J_RANGE[0]}..{J_RANGE[-1]}, N={N}, G={G}\n")
outf.write(f"# gamma  is_real  best_j  best_lsq  best_ev0  log10|ev0|\n")

t0=time.time(); n_done=0
riem_best=[]; fake_best=[]; riem_fail=[]; fake_fp=[]
THRESHOLD=-10  # log|ev0| < -10 → "NEAR-PSD"

print(f"\n{'γ':>12}  {'type':>5}  {'best_j':>7}  {'lsq':>6}  {'best_ev0':>14}  {'log|ev0|':>10}  status")
print("-"*80, flush=True)

with mp.Pool(processes=args.ncores) as pool:
    for res in pool.imap_unordered(worker_jscan, tasks, chunksize=2):
        gamma,is_real,best_j,best_lsq,best_ev=res
        n_done+=1
        log_ev=math.log10(abs(best_ev)) if best_ev not in (float('-inf'),float('nan')) and best_ev!=0 else -99
        typ="REAL" if is_real else "FAKE"
        status=""
        if is_real:
            if log_ev<THRESHOLD: status="OK (near-PSD)"
            else: status="*** MISS ***"; riem_fail.append((gamma,best_ev,best_j))
            riem_best.append((gamma,log_ev))
        else:
            if log_ev<THRESHOLD: status="*** FALSE POS ***"; fake_fp.append((gamma,best_ev,best_j))
            else: status="neg (good)"
            fake_best.append((gamma,log_ev))
        outf.write(f"{gamma:.10f}  {is_real}  {best_j}  {best_lsq}  {best_ev:.6e}  {log_ev:.2f}\n")
        if "***" in status or "near-PSD" in status or n_done%100==0:
            print(f"{gamma:>12.6f}  {typ:>5}  {best_j!s:>7}  {best_lsq!s:>6}  {best_ev:>14.4e}  {log_ev:>10.2f}  {status}", flush=True)
        if n_done%200==0:
            el=time.time()-t0; rate=n_done/el
            print(f"  [{n_done}/{len(tasks)}] {el:.0f}s elapsed, {(len(tasks)-n_done)/rate:.0f}s remaining", flush=True)
            outf.flush()

outf.close()
dt=time.time()-t0

print(f"\n{'='*80}")
print(f"DONE in {dt:.1f}s")
print(f"\nReal zeros: {len(zeros)} tested")
if riem_best:
    log_evs=[x[1] for x in riem_best]
    print(f"  Near-PSD (log|ev0|<-10): {sum(1 for l in log_evs if l<-10)} / {len(zeros)}")
    print(f"  Median log|ev0|: {np.median(log_evs):.2f}")
    print(f"  Best log|ev0|: {min(log_evs):.2f}  (at γ={riem_best[np.argmin(log_evs)][0]:.5f})")
    print(f"  Worst log|ev0|: {max(log_evs):.2f}")
    if riem_fail:
        print(f"  MISSED ({len(riem_fail)} zeros didn't reach -10):")
        for g,e,j in riem_fail[:10]:
            print(f"    γ={g:.5f}: best_j={j}, ev0={e:.3e}")

print(f"\nFakes tested: {len(fakes)}")
if fake_best:
    log_evs=[x[1] for x in fake_best]
    print(f"  False positives (log|ev0|<-10): {len(fake_fp)}")
    print(f"  Median log|ev0|: {np.median(log_evs):.2f}")
    print(f"  Best (least negative) log|ev0|: {min(log_evs):.2f}")
    if fake_fp:
        print(f"  FALSE POSITIVES ({len(fake_fp)}):")
        for g,e,j in fake_fp[:20]:
            print(f"    γ_fake={g:.5f}: best_j={j}, ev0={e:.3e}")

result_summary={
    "J_range":[J_RANGE[0],J_RANGE[-1]],"N":N,"G":G,
    "n_zeros":len(zeros),"n_fakes":len(fakes),
    "zeros_near_psd":sum(1 for g,l in (riem_best or []) if l<-10),
    "fake_false_positives":len(fake_fp),
    "elapsed_s":dt,
    "riem_median_log_ev":float(np.median([x[1] for x in riem_best])) if riem_best else None,
    "fake_median_log_ev":float(np.median([x[1] for x in fake_best])) if fake_best else None,
    "misses":[{"gamma":g,"ev0":e,"best_j":j} for g,e,j in riem_fail[:100]],
    "false_positives":[{"gamma":g,"ev0":e,"best_j":j} for g,e,j in fake_fp[:100]],
}
with open("results/mass_test_jscan.json","w") as f:
    json.dump(result_summary,f,indent=2)
print(f"\nSaved → results/mass_test_jscan.json")

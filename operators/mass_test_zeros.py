#!/usr/bin/env python3
"""
Mass-test: run every LMFDB zero through the CvS operator at j=15, N=64.
Also test 100K fake values (midpoints between consecutive zeros).

Data: ~/projects/riemann/data/zeros_100k.txt (100,000 lines)
Strategy:
  - For each zero γ_k: compute ev0(τ_N) at resonant lsq=round(exp(2π×15/γ_k))
  - For each gap midpoint γ_mid=(γ_k+γ_{k+1})/2: same computation
  - Expected: all real zeros → ev0 ≈ 0 (near-PSD)
             all fakes    → ev0 << 0 (negative)

Multiprocessing: uses all available cores.
Output:
  results/mass_test_zeros.json   (summary + outliers)
  results/mass_test_zeros.out    (full table, streaming)

Usage: python3 operators/mass_test_zeros.py [--nzeros N] [--ncores N]
"""
import math, numpy as np, scipy.linalg, json, sys, time, os
import multiprocessing as mp
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G=60000; N=64; K_BATCH=10000
# Adaptive j: target lsq≈787 for all zeros regardless of γ value
# j = round(γ × log(787) / 2π)  → lsq ≈ round(exp(2π×j/γ)) ≈ 787
L_TARGET=math.log(787)
def adaptive_j(gamma):
    return max(1, round(gamma * L_TARGET / (TWO_PI)))

# Parse args
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--nzeros",type=int,default=100000,help="Number of zeros to test")
parser.add_argument("--ncores",type=int,default=os.cpu_count(),help="Number of worker processes")
parser.add_argument("--batch",type=int,default=100,help="Zeros per batch for progress reporting")
args=parser.parse_args()

N_ZEROS=args.nzeros; N_CORES=args.ncores

# Load zeros
ZEROS_FILE=os.path.expanduser("~/projects/riemann/data/zeros_100k.txt")
print(f"Loading zeros from {ZEROS_FILE}...", flush=True)
with open(ZEROS_FILE) as f:
    all_zeros=[float(line.strip()) for line in f if line.strip()]
zeros=sorted(all_zeros[:N_ZEROS])
print(f"  Loaded {len(zeros)} zeros. Range: {zeros[0]:.4f}..{zeros[-1]:.4f}")

# Build fake values: midpoints
fakes=[(zeros[i]+zeros[i+1])/2 for i in range(len(zeros)-1)]
print(f"  {len(fakes)} fake midpoints")

# Determine lsq range
lsqs_riem=[round(math.exp(TWO_PI*adaptive_j(g)/g)) for g in zeros]
lsqs_fake=[round(math.exp(TWO_PI*adaptive_j(g)/g)) for g in fakes]
max_lsq=max(max(lsqs_riem),max(lsqs_fake))
min_lsq=min(min(lsqs_riem),min(lsqs_fake))
print(f"  lsq range: {min_lsq}..{max_lsq}")

# Sieve primes once
print(f"\nSieving to {max_lsq}...", flush=True)
t0=time.time()
sv=np.ones(max_lsq+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(max_lsq**0.5)+1):
    if sv[i]: sv[i*i::i]=False
primes_arr=np.where(sv)[0]
print(f"  {len(primes_arr)} primes, {time.time()-t0:.1f}s", flush=True)

# Build shared prime data (log, weight) for each UNIQUE lsq
unique_lsqs=sorted(set(lsqs_riem+lsqs_fake))
print(f"  Building prime data for {len(unique_lsqs)} unique lsq values...", flush=True)
t1=time.time()
pd_map={}
for lsq in unique_lsqs:
    ll,ww=[],[]
    for p in primes_arr:
        if p>lsq: break
        lp=math.log(p); pk=p
        while pk<=lsq:
            ll.append(math.log(pk)); ww.append(lp/math.sqrt(pk)); pk*=p
    pd_map[lsq]=(np.array(ll),np.array(ww))
print(f"  Prime data built in {time.time()-t1:.1f}s", flush=True)

def compute_ev0(lsq, N_=None):
    """Compute min eigenvalue of τ_N at given lsq."""
    if N_ is None: N_=N
    L=2*math.log(math.sqrt(lsq)); ns=np.arange(-N_,N_+1); D=len(ns)
    R=TWO_PI*N_/L+50.0; grid=np.linspace(-R,R,G)
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
        S+=np.sin(ph)@wts_p[k0:k1]
        Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd)
    tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])

def worker_fn(args):
    """Worker: returns (gamma, lsq, ev0, is_real)"""
    gamma,lsq,is_real=args
    try:
        ev=compute_ev0(lsq)
        return (gamma,lsq,ev,is_real)
    except Exception as e:
        return (gamma,lsq,float('nan'),is_real)

# Build task list
tasks=[]
for i,(g,lsq) in enumerate(zip(zeros,lsqs_riem)):
    tasks.append((g,lsq,True))
for i,(g,lsq) in enumerate(zip(fakes,lsqs_fake)):
    tasks.append((g,lsq,False))
print(f"\nTotal tasks: {len(tasks)} ({len(zeros)} real + {len(fakes)} fake)")
print(f"Using {N_CORES} worker cores, N={N}, G={G}, adaptive j targeting lsq≈787")
print(f"lsq range for real zeros: {min(lsqs_riem)}..{max(lsqs_riem)}")
print("(Real zeros should have ev0 near 0; fakes should be negative)", flush=True)

os.makedirs("results",exist_ok=True)
outf=open("results/mass_test_zeros.out","w")
outf.write(f"# mass_test_zeros: adaptive_j(lsq~787), N={N}, G={G}\n")
outf.write(f"# gamma  lsq  ev0  is_real  log10|ev0|\n")

t0=time.time()
riem_evs=[]; fake_evs=[]; n_done=0
riem_fails=[]; fake_fails=[]  # zeros that don't test as near-PSD

print(f"\n{'γ':>12}  {'lsq':>8}  {'ev0':>14}  {'log|ev0|':>10}  type  status")
print("-"*70, flush=True)

with mp.Pool(processes=N_CORES) as pool:
    for result in pool.imap_unordered(worker_fn, tasks, chunksize=4):
        gamma,lsq,ev,is_real=result
        n_done+=1
        log_ev=math.log10(abs(ev)) if (ev!=0 and not math.isnan(ev)) else -99
        typ="REAL" if is_real else "FAKE"
        status=""
        if is_real and log_ev < -8: status="OK"
        elif is_real and log_ev >= -8: status="*** FAIL ***"; riem_fails.append((gamma,ev))
        elif not is_real and log_ev > -8: status="*** FALSE POS ***"; fake_fails.append((gamma,ev))
        else: status="OK"
        outf.write(f"{gamma:.10f}  {lsq}  {ev:.6e}  {is_real}  {log_ev:.2f}\n")
        if is_real: riem_evs.append(ev)
        else: fake_evs.append(ev)
        if status not in ("OK","") or n_done % 500 == 0:
            print(f"{gamma:>12.6f}  {lsq:>8}  {ev:>14.4e}  {log_ev:>10.2f}  {typ}  {status}", flush=True)
        if n_done % 1000==0:
            elapsed=time.time()-t0
            rate=n_done/elapsed
            remaining=(len(tasks)-n_done)/rate
            print(f"  [{n_done}/{len(tasks)}] {elapsed:.0f}s elapsed, {remaining:.0f}s remaining", flush=True)
            outf.flush()

outf.close()
dt=time.time()-t0

print(f"\n{'='*70}")
print(f"DONE in {dt:.1f}s ({dt/len(tasks)*1000:.1f}ms/task)")
print(f"\nReal zeros tested: {len(riem_evs)}")
print(f"  Median ev0: {np.median(riem_evs):.4e}")
print(f"  Max ev0: {max(riem_evs):.4e}")
print(f"  Min ev0: {min(riem_evs):.4e}")
print(f"  Near-PSD (log|ev0|<-6): {sum(1 for e in riem_evs if abs(e)<1e-6)} / {len(riem_evs)}")
print(f"\nFake midpoints tested: {len(fake_evs)}")
print(f"  Median ev0: {np.median(fake_evs):.4e}")
print(f"  Max ev0: {max(fake_evs):.4e}")
print(f"  Min ev0: {min(fake_evs):.4e}")
print(f"  False positives (log|ev0|>-6): {len(fake_fails)}")
if fake_fails:
    print(f"  False positive γ values: {[g for g,e in fake_fails[:20]]}")
if riem_fails:
    print(f"\nWARNING: {len(riem_fails)} real zeros FAILED (expected near-PSD but got negative):")
    for g,e in riem_fails[:20]:
        print(f"  γ={g:.6f}: ev0={e:.4e}")

summary={
    "j":J,"N":N,"G":G,
    "n_zeros_tested":len(riem_evs),
    "n_fakes_tested":len(fake_evs),
    "riem_median_ev0":float(np.median(riem_evs)),
    "riem_max_ev0":float(max(riem_evs)),
    "riem_min_ev0":float(min(riem_evs)),
    "riem_near_psd_count":int(sum(1 for e in riem_evs if abs(e)<1e-8)),
    "fake_median_ev0":float(np.median(fake_evs)),
    "fake_false_positive_count":len(fake_fails),
    "riem_fail_count":len(riem_fails),
    "elapsed_s":dt,
    "false_positives":[{"gamma":g,"ev0":e} for g,e in fake_fails[:50]],
    "riem_failures":[{"gamma":g,"ev0":e} for g,e in riem_fails[:50]],
}
with open("results/mass_test_zeros.json","w") as f:
    json.dump(summary,f,indent=2)
print("\nSaved → results/mass_test_zeros.json")
print("      → results/mass_test_zeros.out")

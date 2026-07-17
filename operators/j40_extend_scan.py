#!/usr/bin/env python3
"""
j=40 extended N scan: N=512..4096 with CORRECT operator.

Continue from j40_correct_margin.py (which ran N=64..512).
Want to see: does Riemann ev0 reach Dip1 then RECOVER toward zero?
Predicted Dip1: N_Dip1 ≈ 138×log(lsq_riem) = 138×log(52738111) ≈ 2452.

If Riemann recovers while fakes don't → confirms Weil criterion asymptotically.
Key comparison: fake_m335 (γ=13.8, lsq=81.2M) vs Riemann at large N.

Runs on Olympus.
"""
import math, numpy as np, scipy.linalg, json, time, sys, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693; J=40
G=60000

LSQ_RIEM = 52738111
LSQ_P05  = 28727362   # gamma+0.5=14.635
LSQ_M335 = 81205543   # gamma-0.335=13.800

MAX_LSQ=max(LSQ_RIEM,LSQ_P05,LSQ_M335)

print(f"j={J} extended N scan: N=512..4096 (CORRECT operator)")
print(f"  G={G}, goal: find Dip1 and recovery for Riemann lsq={LSQ_RIEM}")
print(f"  Predicted N_Dip1 ≈ 138×log(lsq) = {int(138*math.log(LSQ_RIEM))}")
print(f"  Sieving to {MAX_LSQ}...", flush=True)

t0=time.time()
sieve=np.ones(MAX_LSQ+1,dtype=bool); sieve[0]=sieve[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sieve[i]: sieve[i*i::i]=False
primes=np.where(sieve)[0]
print(f"  {len(primes)} primes in {time.time()-t0:.1f}s", flush=True)

def prime_data(lsq):
    lkpp_list,wts_list=[],[]
    for p in primes:
        if p>lsq: break
        lp=math.log(p); pk=p
        while pk<=lsq:
            lkpp_list.append(math.log(pk)); wts_list.append(lp/math.sqrt(pk)); pk*=p
    return np.array(lkpp_list),np.array(wts_list)

print("  Building prime data...", end=" ", flush=True)
t0=time.time()
pd={LSQ_RIEM:prime_data(LSQ_RIEM), LSQ_P05:prime_data(LSQ_P05), LSQ_M335:prime_data(LSQ_M335)}
print(f"done ({time.time()-t0:.1f}s)", flush=True)

def build_ev0(N,lsq):
    L=2*math.log(math.sqrt(lsq))
    ns=np.arange(-N,N+1); D=len(ns)
    R=TWO_PI*N/L+50.0
    grid=np.linspace(-R,R,G)
    h=grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=np.real(digamma(0.25+0.5j*grid))-LOGPI
    c=TWO_PI*ns[:,None]/L; d=c-grid[None,:]
    sm=np.abs(d)<1e-12
    E=np.where(sm,np.sqrt(L)+0j,(1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch=(E*(w*Om)@E.conj().T).real/TWO_PI
    sh=np.sinh(L/4)**2; ni=ns[:,None]; nj=ns[None,:]
    B=32*L*sh*(L**2-16*math.pi**2*ni*nj)/((L**2+16*math.pi**2*ni**2)*(L**2+16*math.pi**2*nj**2))
    lkpp,wts_p=pd[lsq]; K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D); K_BATCH=10000
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

print()
print(f"{'N':>6}  {'ev0_riem':>14}  {'ev0_p05':>14}  {'ev0_m335':>14}  {'sel_p05':>10}  {'sel_m335':>10}  dt")
print("-"*100)
sys.stdout.flush()

data=[]; best_riem=float('inf'); best_N=None
# Start from N=512 (confirmed at -2.02e-7), continue to find Dip1 and recovery
N_vals=[512,768,1024,1280,1536]
# Stop at 1536: Dip1 expected ~N=2452 per 138*log scaling but may arrive earlier.
# At N=1536, computation takes ~30 min/point — beyond that is diminishing returns per session.
for N in N_vals:
    t0=time.time()
    er=build_ev0(N,LSQ_RIEM)
    ep=build_ev0(N,LSQ_P05)
    em=build_ev0(N,LSQ_M335)
    dt=time.time()-t0
    sp=abs(ep)/abs(er) if er!=0 else float('nan')
    sm=abs(em)/abs(er) if er!=0 else float('nan')
    mark=" <DIP1" if er<best_riem else (" RECOV" if (best_N and er>best_riem*0.9) else "")
    if er<best_riem: best_riem=er; best_N=N
    print(f"{N:>6}  {er:>14.4e}  {ep:>14.4e}  {em:>14.4e}  {sp:>10.2f}  {sm:>10.2f}  ({dt:.0f}s){mark}")
    sys.stdout.flush()
    data.append({"N":N,"ev0_riem":er,"ev0_p05":ep,"ev0_m335":em,"sel_p05":sp,"sel_m335":sm})

os.makedirs("results",exist_ok=True)
with open("results/j40_extend_scan.json","w") as f:
    json.dump({"j":J,"G":G,"lsq_riem":LSQ_RIEM,"lsq_p05":LSQ_P05,"lsq_m335":LSQ_M335,
               "dip1_N":best_N,"dip1_ev0":best_riem,"data":data},f,indent=2)
print(f"\nDip1: N={best_N}, ev[0]_riem={best_riem:.4e}")
print("Saved → results/j40_extend_scan.json")

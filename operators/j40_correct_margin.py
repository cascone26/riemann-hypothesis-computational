#!/usr/bin/env python3
"""
j=40 CORRECT operator selectivity — small N regime.

KEY FINDING (j=20 analysis): selectivity lives at SMALL N (N≤256).
At large N, arch dominates and all lsq converge to same ev[0].
The 5782× ratio must be at small N with j=40's massive prime sum (~3.1M primes).

Strategy: N=64,128,256,512 only, G=60000 for accuracy.
Compare Riemann lsq vs fake lsq at each N.

Primes precomputed ONCE (sieve to max_lsq=81M, then filtered per lsq).
Runs on Olympus (cobo@100.124.234.43).
"""
import math, numpy as np, scipy.linalg, json, time, sys, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693; J=40
G=60000

# j=40 fake gammas
# gamma1+0.5=14.635: lsq=round(exp(2pi*40/14.635))
# gamma1-0.335=13.800: lsq=round(exp(2pi*40/13.800)) — the "13.8" fake
# gamma1+1.0=15.135: lsq=round(exp(2pi*40/15.135))
LSQ_RIEM = round(math.exp(TWO_PI*J/GAMMA1))
LSQ_P05  = round(math.exp(TWO_PI*J/(GAMMA1+0.5)))
LSQ_M335 = round(math.exp(TWO_PI*J/(GAMMA1-0.335)))  # "13.8" fake
LSQ_P10  = round(math.exp(TWO_PI*J/(GAMMA1+1.0)))

MAX_LSQ = max(LSQ_RIEM, LSQ_P05, LSQ_M335, LSQ_P10)

print(f"j={J} CORRECT operator selectivity — SMALL N regime")
print(f"  G={G} (large for accuracy)")
print(f"  lsq_riem  = {LSQ_RIEM}")
print(f"  lsq_p05   = {LSQ_P05}   (gamma1+0.5={GAMMA1+0.5:.3f})")
print(f"  lsq_m335  = {LSQ_M335}  (gamma1-0.335=13.800, '13.8 fake')")
print(f"  lsq_p10   = {LSQ_P10}   (gamma1+1.0={GAMMA1+1.0:.3f})")
print(f"  Sieving to {MAX_LSQ}...")
sys.stdout.flush()

t0=time.time()
sieve=np.ones(MAX_LSQ+1,dtype=bool); sieve[0]=sieve[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sieve[i]: sieve[i*i::i]=False
primes=np.where(sieve)[0]
print(f"  {len(primes)} primes found in {time.time()-t0:.1f}s")
sys.stdout.flush()

# Precompute (lkpp, wts_p) for each lsq
def prime_data(lsq):
    lkpp_list, wts_list = [], []
    for p in primes:
        if p > lsq: break
        lp = math.log(p); pk = p
        while pk <= lsq:
            lkpp_list.append(math.log(pk)); wts_list.append(lp/math.sqrt(pk)); pk*=p
    return np.array(lkpp_list), np.array(wts_list)

print("  Building prime data for each lsq...")
t0=time.time()
pd = {
    LSQ_RIEM: prime_data(LSQ_RIEM),
    LSQ_P05:  prime_data(LSQ_P05),
    LSQ_M335: prime_data(LSQ_M335),
    LSQ_P10:  prime_data(LSQ_P10),
}
print(f"  Done ({time.time()-t0:.1f}s). Prime terms: "
      f"riem={len(pd[LSQ_RIEM][0])}, p05={len(pd[LSQ_P05][0])}, "
      f"m335={len(pd[LSQ_M335][0])}, p10={len(pd[LSQ_P10][0])}")
sys.stdout.flush()

def build_ev0(N, lsq):
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
    lkpp,wts_p=pd[lsq]
    K=len(lkpp); S=np.zeros(D); Sd=np.zeros(D)
    K_BATCH=10000
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
print(f"{'N':>6}  {'ev0_riem':>14}  {'ev0_p05':>14}  {'ev0_m335':>14}  {'ev0_p10':>14}  "
      f"{'sel_p05':>10}  {'sel_m335':>10}  {'sel_p10':>10}  dt")
print("-"*125)
sys.stdout.flush()

data=[]; N_vals=[64,128,256,512]
for N in N_vals:
    t0=time.time()
    er=build_ev0(N,LSQ_RIEM)
    ep=build_ev0(N,LSQ_P05)
    em=build_ev0(N,LSQ_M335)
    ep10=build_ev0(N,LSQ_P10)
    dt=time.time()-t0
    # sel = |ev0_fake|/|ev0_riem|, > 1 means riem closer to 0 (better Weil positivity)
    sp  = abs(ep)/abs(er)   if er!=0 else float('nan')
    sm  = abs(em)/abs(er)   if er!=0 else float('nan')
    sp10= abs(ep10)/abs(er) if er!=0 else float('nan')
    print(f"{N:>6}  {er:>14.4e}  {ep:>14.4e}  {em:>14.4e}  {ep10:>14.4e}  "
          f"{sp:>10.2f}  {sm:>10.2f}  {sp10:>10.2f}  ({dt:.0f}s)")
    sys.stdout.flush()
    data.append({"N":N,"ev0_riem":er,"ev0_p05":ep,"ev0_m335":em,"ev0_p10":ep10,
                 "sel_p05":sp,"sel_m335":sm,"sel_p10":sp10})

print()
print("=== Summary: does j=40 show 5782× selectivity? ===")
print("  sel = |ev0_fake|/|ev0_riem|  — >1 means Riemann closer to zero (CORRECT Weil behavior)")
for d in data:
    print(f"  N={d['N']:>4}: sel_m335={d['sel_m335']:>10.2f}  (vs claimed 5782× from abc_credence)")

os.makedirs("results",exist_ok=True)
with open("results/j40_correct_margin.json","w") as f:
    json.dump({"j":J,"G":G,"lsq_riem":LSQ_RIEM,"lsq_p05":LSQ_P05,
               "lsq_m335":LSQ_M335,"lsq_p10":LSQ_P10,"data":data},f,indent=2)
print("\nSaved → results/j40_correct_margin.json")

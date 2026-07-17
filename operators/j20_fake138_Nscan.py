#!/usr/bin/env python3
"""
j=20 N-sweep: verify that fake_13.8 false-negative at N=64 is a finite-N artifact.

At N=64, j=20: ev0_riem=-9.50e-10, ev0_fake=-4.23e-12, sel=0.004 (false negative).
Hypothesis: at larger N, riem recovers toward zero faster than fake → sel→∞.

Mechanism: fake has lsq=9010 > riem lsq=7262 → more primes → larger P → better cancellation at
small N. But Weil criterion demands ONLY Riemann lsq shows true near-PSD in the N→∞ limit.

Tests: N=64,128,256,384,512, G=30000, both lsqs pre-sieved.
"""
import math, numpy as np, scipy.linalg, time, sys
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693
GAMMA_FAKE=13.800; J=20; G=30000

LSQ_RIEM=round(math.exp(TWO_PI*J/GAMMA1))  # 7262
LSQ_FAKE=round(math.exp(TWO_PI*J/GAMMA_FAKE))  # 9010

print(f"j={J} N-sweep: riem lsq={LSQ_RIEM}, fake13.8 lsq={LSQ_FAKE}")
print(f"Goal: show sel_fake13.8 grows with N (false negative at N=64 is finite-N artifact)\n")

# Precompute prime data for both lsqs
def sieve_lsq(lsq):
    sv=np.ones(lsq+1,dtype=bool); sv[0]=sv[1]=False
    for i in range(2,int(lsq**0.5)+1):
        if sv[i]: sv[i*i::i]=False
    lkpp,wts=[],[]
    for p in range(2,lsq+1):
        if not sv[p]: continue
        lp=math.log(p); pk=p
        while pk<=lsq:
            lkpp.append(math.log(pk)); wts.append(lp/math.sqrt(pk)); pk*=p
    return np.array(lkpp), np.array(wts)

print("Sieving primes...", end=" ", flush=True)
pd={LSQ_RIEM: sieve_lsq(LSQ_RIEM), LSQ_FAKE: sieve_lsq(LSQ_FAKE)}
print(f"done. riem K={len(pd[LSQ_RIEM][0])}, fake K={len(pd[LSQ_FAKE][0])}")

def build_ev0(N, lsq):
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
    lkpp,wts_p=pd[lsq]; K=len(lkpp)
    S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,10000):
        k1=min(k0+10000,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]
        Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd)
    tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])

print(f"\n{'N':>5}  {'ev0_riem':>14}  {'ev0_fake':>14}  {'sel':>12}  dt")
print("-"*65)
sys.stdout.flush()

for N in [64, 128, 192, 256, 384, 512]:
    t0=time.time()
    er=build_ev0(N, LSQ_RIEM)
    ef=build_ev0(N, LSQ_FAKE)
    dt=time.time()-t0
    sel=abs(ef)/abs(er) if er!=0 and ef!=0 else float('nan')
    note=""
    if sel < 1: note=" ← FALSE NEG"
    elif sel > 100: note=" ← GOOD"
    print(f"{N:>5}  {er:>14.4e}  {ef:>14.4e}  {sel:>12.2f}  ({dt:.0f}s){note}")
    sys.stdout.flush()

print("\nDone.")

#!/usr/bin/env python3
"""Single-point j=40, N=1536 — continuation of j40_extend_scan."""
import math, numpy as np, scipy.linalg, time, sys
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693; G=60000; N=1536
LSQ_RIEM=52738111; LSQ_P05=28727362; LSQ_M335=81205543
MAX_LSQ=max(LSQ_RIEM,LSQ_P05,LSQ_M335)

print(f"j=40, N={N}: single point, G={G}")
t0=time.time()
sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sv[i]: sv[i*i::i]=False
primes=np.where(sv)[0]
print(f"  {len(primes)} primes in {time.time()-t0:.1f}s", flush=True)

def prime_data(lsq):
    ll,ww=[],[]
    for p in primes:
        if p>lsq: break
        lp=math.log(p); pk=p
        while pk<=lsq:
            ll.append(math.log(pk)); ww.append(lp/math.sqrt(pk)); pk*=p
    return np.array(ll),np.array(ww)

print("  Building prime data...", flush=True)
pd={lsq: prime_data(lsq) for lsq in [LSQ_RIEM,LSQ_P05,LSQ_M335]}

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
    lkpp,wts_p=pd[lsq]; K=len(lkpp); S=np.zeros(D); Sd=np.zeros(D)
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

print(f"\n  N={N}: computing...", flush=True)
t0=time.time()
er=build_ev0(N,LSQ_RIEM); print(f"  riem: {er:.4e}", flush=True)
ep=build_ev0(N,LSQ_P05); print(f"  p05:  {ep:.4e}", flush=True)
em=build_ev0(N,LSQ_M335); print(f"  m335: {em:.4e}", flush=True)
dt=time.time()-t0
sp=abs(ep)/abs(er); sm=abs(em)/abs(er)
print(f"\n  N={N}: ev0_riem={er:.4e}, ev0_p05={ep:.4e}, ev0_m335={em:.4e}")
print(f"  sel_p05={sp:.2f}, sel_m335={sm:.2f}, dt={dt:.0f}s")

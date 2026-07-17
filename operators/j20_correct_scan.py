#!/usr/bin/env python3
"""
j=20 N-scan with CORRECT operator (Sj-Si sign fixed).

At N=64: Riemann ev[0]=-9.5e-10 (nearly PSD), fake+0.5 ev[0]=-1.7e-5 (negative).
At N=512: Riemann ev[0]=-2.54e-3 (deeper), fake+0.5=-7.05e-4 (recovering).

Find: where is the Riemann Dip1? Does it recover toward 0 after?
Compare with fake lsq at each N.
"""
import math, numpy as np, scipy.linalg, json, time, sys, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
GAMMA1=14.134725141734693
J=20
G=3000; K_BATCH=3000

LSQ_RIEM = round(math.exp(TWO_PI*J/GAMMA1))   # 7262
LSQ_P05  = round(math.exp(TWO_PI*J/(GAMMA1+0.5)))  # ~5305
LSQ_M05  = round(math.exp(TWO_PI*J/(GAMMA1-0.5)))  # ~9121

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
    sieve=np.ones(lsq+1,dtype=bool); sieve[0]=sieve[1]=False
    for i in range(2,int(lsq**0.5)+1):
        if sieve[i]: sieve[i*i::i]=False
    lkpp,wts_p=[],[]
    for p in range(2,lsq+1):
        if not sieve[p]: continue
        pk=p
        while pk<=lsq:
            lkpp.append(math.log(pk)); wts_p.append(math.log(p)/math.sqrt(pk)); pk*=p
    lkpp=np.array(lkpp); wts_p=np.array(wts_p)
    K=len(lkpp); S=np.zeros(D); Sd=np.zeros(D)
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

print(f"j={J}: CORRECT operator N-scan")
print(f"  lsq_riem={LSQ_RIEM}  lsq_p05={LSQ_P05}  lsq_m05={LSQ_M05}")
print()
print(f"{'N':>6}  {'ev0_riem':>12}  {'ev0_p05':>12}  {'ev0_m05':>12}  {'R_p05':>10}  {'R_m05':>10}  dt")
print("-"*90)
sys.stdout.flush()

data=[]
N_vals=[64,128,256,512,768,1024,1280,1536,1792,2048,2560,3072,4096]
for N in N_vals:
    t0=time.time()
    er=build_ev0(N,LSQ_RIEM)
    ep=build_ev0(N,LSQ_P05)
    em=build_ev0(N,LSQ_M05)
    dt=time.time()-t0
    rp=abs(er)/abs(ep) if ep!=0 else float('nan')
    rm=abs(er)/abs(em) if em!=0 else float('nan')
    print(f"{N:>6}  {er:>12.4e}  {ep:>12.4e}  {em:>12.4e}  {rp:>10.4f}  {rm:>10.4f}  ({dt:.1f}s)")
    sys.stdout.flush()
    data.append({"N":N,"ev0_riem":er,"ev0_p05":ep,"ev0_m05":em,"R_p05":rp,"R_m05":rm})

os.makedirs("results",exist_ok=True)
with open("results/j20_correct_scan.json","w") as f:
    json.dump({"j":J,"lsq_riem":LSQ_RIEM,"lsq_p05":LSQ_P05,"lsq_m05":LSQ_M05,"data":data},f,indent=2)
print("\nSaved → results/j20_correct_scan.json")

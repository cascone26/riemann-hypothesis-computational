#!/usr/bin/env python3
"""
high_N_confirm.py — Extend N to 512 for key points.

lsq=787: γ₁ resonance, confirmed STABLE at N=256. Does it approach 0?
lsq=797: borderline STABLE at N=256 (dip at N=192). Resolve: still stable or diverges?
lsq=967: γ₂ resonance, STABLE at N=128. Confirm at N=384.
lsq=926: known DIVERGING control. Confirm.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

TEST_LSQS=[
    (787, "γ₁ resonance", [256,320,384,448,512]),
    (797, "borderline stable?", [256,320,384,448,512]),
    (967, "γ₂ resonance", [192,256,320,384]),
    (926, "known diverging control", [128,192,256]),
]

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
    p_in=primes[primes<=lsq]; ll=[]; ww=[]
    for p in p_in:
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

results={}
for lsq,label,N_vals in TEST_LSQS:
    print(f"\n=== lsq={lsq} ({label}) ===")
    row=[]
    for N in N_vals:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
        marker=" STABLE" if log_e<-10 else " DIVERGED"
        print(f"  N={N:>3}: ev0={e:>14.6e}  log={log_e:>7.2f}  {dt:.1f}s{marker}", flush=True)
    results[str(lsq)]={"label":label,"data":row}

print("\n\n=== SUMMARY ===")
print(f"{'lsq':>5}  {'label':<28}  N-values and log|ev0|")
for lsq,label,_ in TEST_LSQS:
    rd={r["N"]:r["log10_abs_ev0"] for r in results[str(lsq)]["data"]}
    vals=" ".join(f"N{N}:{rd[N]:+.1f}" for N in sorted(rd))
    print(f"{lsq:>5}  {label:<28}  {vals}")

os.makedirs("results",exist_ok=True)
with open("results/high_N_confirm.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/high_N_confirm.json")

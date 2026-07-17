#!/usr/bin/env python3
"""light_checks.py — Smaller N tests (no N>256 to avoid memory kills)."""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900
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
    S=np.zeros(D); Sd=np.zeros(D)
    for k0 in range(0,K,10000):
        k1=min(k0+10000,K)
        ph=TWO_PI*np.outer(ns,lkpp[k0:k1])/L
        S+=np.sin(ph)@wts_p[k0:k1]; Sd+=np.cos(ph)@(2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    dn=ns[:,None]-ns[None,:]; Si=S[:,None]; Sj=S[None,:]
    with np.errstate(divide='ignore',invalid='ignore'):
        P=np.where(dn!=0,(Sj-Si)/(math.pi*dn),0.0)
    np.fill_diagonal(P,Sd); tau=(arch+B-P); tau=(tau+tau.T)/2
    return float(scipy.linalg.eigvalsh(tau,subset_by_index=[0,0])[0])

TESTS=[
    ("lsq=967 N=256,320",  967, [256, 320]),
    ("γ₈ j=43 lsq=508 ext", 508, [192, 256]),
    ("γ₅ j=32 lsq=450",  450, [32,64,96,128,192]),
    ("γ₅ j=38 lsq=1407", 1407,[32,64,128,192]),  # high j for γ₅
]

results={}
for label,lsq,N_vals in TESTS:
    print(f"\n=== {label} ===")
    row=[]
    for N in N_vals:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"log10_abs_ev0":log_e})
        marker=" STABLE" if log_e<-10 else ""
        print(f"  N={N:>3}: log={log_e:>7.2f}  {dt:.1f}s{marker}", flush=True)
    results[label]=row

print("\n\n=== SUMMARY ===")
for label,lsq,N_vals in TESTS:
    row=results[label]
    rd={r["N"]:r["log10_abs_ev0"] for r in row}
    nmin,nmax=min(rd),max(rd); delta=rd[nmax]-rd[nmin]
    all_stable=all(v<-10 for v in rd.values())
    verdict=("STABLE" if abs(delta)<2 and all_stable else ("PSEUDO" if abs(delta)<2 else "DIVERGING"))
    vals=" ".join(f"N{N}:{rd[N]:.1f}" for N in sorted(rd))
    print(f"  {label}: {vals} → {verdict}")

os.makedirs("results",exist_ok=True)
with open("results/light_checks.json","w") as f: json.dump(results,f,indent=2)
print("Saved → results/light_checks.json")

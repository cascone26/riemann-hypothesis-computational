#!/usr/bin/env python3
"""
gamma4_jscan.py — Find stable resonance for γ₄=30.424876.

j=32 → lsq=741 showed DIVERGING at N=128. Test other j values.
Candidates: j=31→603, j=33→911, j=34→1120

If γ₄ has a stable resonance at any j, it confirms universality.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900
GAMMA4=30.424876125859513

# Compute candidates
candidates=[]
for j in range(30, 37):
    lsq=round(math.exp(TWO_PI*j/GAMMA4))
    candidates.append((j,lsq))

print("γ₄=30.424876 j-scan:")
for j,lsq in candidates:
    print(f"  j={j}: lsq={lsq}")
print()

MAX_LSQ=1200
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

N_VALUES=[32,64,96,128,192]
results={}

for j,lsq in candidates:
    if lsq>MAX_LSQ:
        print(f"Skipping j={j} lsq={lsq} (>MAX_LSQ={MAX_LSQ})")
        continue
    label=f"γ₄ j={j}"
    print(f"\n=== lsq={lsq} ({label}) ===")
    row=[]
    for N in N_VALUES:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
        marker=" *** DEEP ***" if log_e<-10 else ""
        print(f"  N={N:>3}: ev0={e:>12.4e}  log={log_e:>7.2f}  {dt:.1f}s{marker}", flush=True)
    results[str(lsq)]={"j":j,"gamma4":GAMMA4,"data":row}

print("\n\n=== SUMMARY: γ₄ j-scan ===")
print(f"{'j':>4}  {'lsq':>5}  {'N=32':>7}  {'N=64':>7}  {'N=128':>7}  {'N=192':>7}  verdict")
print("-"*65)
for j,lsq in candidates:
    if str(lsq) not in results: continue
    rd={r["N"]:r["log10_abs_ev0"] for r in results[str(lsq)]["data"]}
    v32=rd.get(32,99); v64=rd.get(64,99); v128=rd.get(128,99); v192=rd.get(192,99)
    delta=v192-v32 if 192 in rd else v128-v32
    verdict="STABLE" if abs(delta)<2 else "DIVERGING"
    known=" (known)" if j==32 else ""
    print(f"{j:>4}  {lsq:>5}  {v32:>7.2f}  {v64:>7.2f}  {v128:>7.2f}  {v192:>7.2f}  {verdict}{known}")

os.makedirs("results",exist_ok=True)
with open("results/gamma4_jscan.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/gamma4_jscan.json")

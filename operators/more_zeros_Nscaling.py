#!/usr/bin/env python3
"""
more_zeros_Nscaling.py — Test N-scaling for γ₅..γ₁₀.

For each zero, test the j value that gives lsq nearest ~800-900.
Key question: is N-scaling stability universal across Riemann zeros?

Known: γ₁(j=15,lsq=787), γ₂(j=23,lsq=967), γ₃(j=27,lsq=883) all STABLE.
       γ₄(j=32,lsq=741) DIVERGING (anomaly — other j tested in gamma4_jscan.py).

Also tests: γ₅ at both j=35(lsq=794) [already in cluster data, control]
            and j=36(lsq=961) [new].
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

ZEROS=[
    (5,  32.935061587739189, [(35,794),(36,961)]),
    (6,  37.586178158825671, [(40,802),(41,948)]),
    (7,  40.918719012147495, [(44,859),(45,1002)]),
    (8,  43.327073280914999, [(46,789),(47,912)]),
    (9,  48.005150881167159, [(51,792),(52,903)]),
    (10, 49.773832477672302, [(53,805),(54,913)]),
]

MAX_LSQ=1100
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

N_VALUES=[32,64,96,128]
results={}

for k,gamma,js in ZEROS:
    for j,lsq in js:
        if lsq>MAX_LSQ:
            print(f"Skipping γ_{k} j={j} lsq={lsq} (>MAX_LSQ)")
            continue
        key=f"g{k}_j{j}_lsq{lsq}"
        label=f"γ_{k}={gamma:.3f} j={j}"
        print(f"\n=== {label} → lsq={lsq} ===")
        row=[]
        for N in N_VALUES:
            t0=time.time()
            e=compute_ev0(lsq,N)
            dt=time.time()-t0
            log_e=math.log10(abs(e)) if abs(e)>0 else -99
            row.append({"N":N,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
            marker=" NEAR-PSD" if log_e<-10 else ""
            print(f"  N={N:>3}: ev0={e:>12.4e}  log={log_e:>7.2f}  {dt:.1f}s{marker}", flush=True)
        results[key]={"k":k,"gamma":gamma,"j":j,"lsq":lsq,"data":row}

print("\n\n=== SUMMARY: γ₅..γ₁₀ N-scaling ===")
print(f"{'zero':>8}  {'j':>4}  {'lsq':>5}  {'N=32':>7}  {'N=64':>7}  {'N=128':>7}  Δlog  verdict")
print("-"*75)
for k,gamma,js in ZEROS:
    for j,lsq in js:
        key=f"g{k}_j{j}_lsq{lsq}"
        if key not in results: continue
        rd={r["N"]:r["log10_abs_ev0"] for r in results[key]["data"]}
        v32=rd.get(32,99); v64=rd.get(64,99); v128=rd.get(128,99)
        delta=v128-v32
        verdict="STABLE" if abs(delta)<2 else "DIVERGING"
        star=" ***" if verdict=="STABLE" else ""
        print(f"  γ_{k}({gamma:.3f})  j={j:>2}  {lsq:>5}  {v32:>7.2f}  {v64:>7.2f}  {v128:>7.2f}  {delta:>+5.2f}  {verdict}{star}")

os.makedirs("results",exist_ok=True)
with open("results/more_zeros_Nscaling.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/more_zeros_Nscaling.json")

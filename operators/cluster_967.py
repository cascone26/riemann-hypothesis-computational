#!/usr/bin/env python3
"""
cluster_967.py — N-scaling cluster collapse for γ₂ neighborhood.

γ₂=21.022, j=23 → lsq=967. Does the same cluster collapse happen here?
cluster_N_scaling.py showed lsq=787 uniquely stable in [785-799].
This tests whether lsq=967 is uniquely stable in its neighborhood.

Also: γ₅ at j=36 → lsq=961. Does lsq=961 appear as a secondary stability?
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

CLUSTER_LSQS=list(range(960, 978))  # 960..977
N_VALUES=[32, 64, 96, 128, 192, 256]

MAX_LSQ=max(CLUSTER_LSQS)+5
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

print("CLUSTER N-SCALING: γ₂ neighborhood [960-977]")
print(f"γ₂=21.022, j=23 → lsq=967  (resonance under test)")
print(f"γ₅=32.935, j=36 → lsq=961  (secondary candidate)")
print(f"N values: {N_VALUES}")
print()

results={}
for lsq in CLUSTER_LSQS:
    print(f"\n=== lsq={lsq} ===")
    row=[]
    for N in N_VALUES:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
        marker=" ***" if log_e<-10 else ""
        print(f"  N={N:>3}: log|ev0|={log_e:>7.2f}  {dt:.1f}s{marker}", flush=True)
    results[lsq]=row

GAMMA2=21.022039638771554; J2=23
res2=round(math.exp(TWO_PI*J2/GAMMA2))
GAMMA5=32.935061587739189; J5=36
res5=round(math.exp(TWO_PI*J5/GAMMA5))
print(f"\n\n=== SUMMARY: log₁₀|ev₀| ===")
print(f"(γ₂ resonant lsq at j={J2}: {res2})")
print(f"(γ₅ resonant lsq at j={J5}: {res5})")
print()
header=f"{'lsq':>5}"
for N in N_VALUES: header+=f"  N={N:>3}"
header+="  verdict  note"
print(header)
print("-"*(len(header)+15))

for lsq in CLUSTER_LSQS:
    rd={r["N"]:r["log10_abs_ev0"] for r in results[lsq]}
    v32=rd.get(32,99); v256=rd.get(256,99); v128=rd.get(128,99)
    delta=v256-v32 if 256 in rd else v128-v32
    verdict="STABLE" if abs(delta)<2 else "DIVERGING"
    note=""
    if lsq==res2: note=f" <-- γ₂ RESONANCE"
    if lsq==res5: note=f" <-- γ₅ j={J5} candidate"
    if abs(lsq-res5)<=1 and lsq!=res5: note=f" <-- near γ₅ resonance"
    line=f"{lsq:>5}"
    for N in N_VALUES: line+=f"  {rd.get(N,99):>5.2f}"
    line+=f"  {verdict}{note}"
    print(line)

os.makedirs("results",exist_ok=True)
with open("results/cluster_967.json","w") as f:
    json.dump({"N_values":N_VALUES,"results":{str(lsq):v for lsq,v in results.items()}},f,indent=2)
print(f"\nSaved → results/cluster_967.json")

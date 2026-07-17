#!/usr/bin/env python3
"""
cluster_N_scaling.py — N-scaling for all lsq in the γ₁ cluster.

Key test: within the near-PSD cluster [785-799] at N=64, which lsq values
remain near-PSD as N → 256? If only the Riemann resonant lsq=787-789 stays
stable, this is a genuine zero-discriminator operating via N-scaling.

Tests lsq = 785, 786, 787, 788, 789, 790, 795, 797, 798, 799
at N = 32, 64, 96, 128, 192, 256

γ₁ = 14.1347: resonant lsq ≈ 787-789 at j=15
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

# All lsq in the γ₁ cluster region
CLUSTER_LSQS=[
    785, 786, 787, 788, 789, 790,  # dense cluster (all near-PSD at N=64)
    791, 792, 793, 794,             # gap region (NOT near-PSD at N=64)
    795, 796, 797, 798, 799,        # next cluster (near-PSD at N=64)
]
N_VALUES=[32, 64, 96, 128, 192, 256]

MAX_LSQ=max(CLUSTER_LSQS)+10
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
    ll=[]; ww=[]
    for p in primes:
        if p>lsq: break
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

print("CLUSTER N-SCALING: γ₁ cluster [785-799]")
print(f"N values: {N_VALUES}")
print(f"γ₁=14.1347 resonance: lsq ≈ 787-789 (j=15)")
print()

results={}
for lsq in CLUSTER_LSQS:
    print(f"\n=== lsq={lsq} ===")
    print(f"{'N':>5}  {'ev0':>14}  {'log|ev0|':>10}  {'time':>7}")
    print("-"*40)
    row=[]
    for N in N_VALUES:
        G=G_PER_N*N
        t0=time.time()
        e=compute_ev0(lsq,N,G)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
        marker=""
        if log_e>-10: marker=" *** NEGATIVE ***"
        print(f"{N:>5}  {e:>14.6e}  {log_e:>10.2f}  {dt:>6.1f}s{marker}", flush=True)
    results[lsq]=row

print("\n\n=== SUMMARY: log₁₀|ev0| ===")
header=f"{'lsq':>5}"
for N in N_VALUES: header+=f"  N={N:>3}"
header+="  verdict"
print(header)
print("-"*(len(header)+10))

GAMMA1=14.1347251
J=15
res_lsq=round(math.exp(TWO_PI*J/GAMMA1))
print(f"  (γ₁ resonant lsq at j={J}: {res_lsq})")
print()

for lsq in CLUSTER_LSQS:
    row_data={r["N"]:r["log10_abs_ev0"] for r in results[lsq]}
    v32=row_data.get(32,99); v64=row_data.get(64,99); v128=row_data.get(128,99); v256=row_data.get(256,99)
    delta=v256-v32 if 256 in row_data else v128-v32
    verdict="STABLE" if abs(delta)<2 else "DIVERGING"
    resonant=" <-- γ₁ RESONANCE" if abs(lsq-res_lsq)<=2 else ""
    line=f"{lsq:>5}"
    for N in N_VALUES: line+=f"  {row_data.get(N,99):>5.2f}"
    line+=f"  {verdict}{resonant}"
    print(line)

os.makedirs("results",exist_ok=True)
with open("results/cluster_N_scaling.json","w") as f:
    json.dump({"N_values":N_VALUES,"results":{str(lsq):v for lsq,v in results.items()}},f,indent=2)
print(f"\nSaved → results/cluster_N_scaling.json")

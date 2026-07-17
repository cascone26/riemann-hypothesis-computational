#!/usr/bin/env python3
"""
N_scaling_v2.py — Extended N-scaling study.

Key finding from v1: lsq=787 and lsq=967 (Riemann zero resonances) MAINTAIN
ev0 ~ -10^{-11} up to N=128, while lsq=926, 779, 500 DIVERGE to -10^{-4}.

This version:
1. Extends N to 256 for key points (787 vs 926 only)
2. Adds γ₃ resonance (lsq=883), γ₄ resonance (lsq=741)
3. Adds neighboring non-resonant lsq values (788,789 near γ₁; 779 already done)
4. Tests a second "fake" (lsq_j=15 for γ=14.5, lsq=round(exp(2π×15/14.5))=680)

Goal: confirm that Riemann-zero resonances are uniquely stable under N-scaling.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

TEST_LSQS=[
    (787, "γ₁ resonance (j=15, γ₁=14.1347)"),
    (967, "γ₂ resonance (j=23, γ₂=21.022)"),
    (883, "γ₃ resonance (j=27, γ₃=25.011)"),
    (741, "γ₄ resonance (j=32, γ₄=30.425)"),
    (788, "γ₁ neighbor +1"),
    (779, "gap: lsq=779 (between clusters)"),
    (926, "fake: j=15, γ=13.8 → lsq=926"),
    (680, "fake: j=15, γ=14.5 → lsq=680"),
]
N_VALUES_STANDARD=[8, 16, 32, 64, 96, 128]
N_VALUES_EXTENDED=[8, 16, 32, 64, 96, 128, 192, 256]

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

print("N-SCALING STUDY v2: Extended test")
print(f"Grid = {G_PER_N} × N")
print()

results={}
for lsq,label in TEST_LSQS:
    # Use extended range only for key comparison points
    use_extended=(lsq in [787, 926])
    N_vals=N_VALUES_EXTENDED if use_extended else N_VALUES_STANDARD
    print(f"\n=== lsq={lsq} ({label}) ===")
    print(f"{'N':>5}  {'G':>7}  {'ev0':>14}  {'log|ev0|':>10}  {'time':>7}")
    print("-"*50)
    row=[]
    for N in N_vals:
        G=G_PER_N*N
        t0=time.time()
        e=compute_ev0(lsq,N,G)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"G":G,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
        print(f"{N:>5}  {G:>7}  {e:>14.6e}  {log_e:>10.2f}  {dt:>6.1f}s", flush=True)
    results[lsq]={"label":label,"data":row}

print("\n\n=== SUMMARY TABLE ===")
print(f"{'lsq':>6}  {'label':<45}  N=32   N=64   N=128  N=256")
print("-"*100)
for lsq,label in TEST_LSQS:
    row_data={r["N"]:r["log10_abs_ev0"] for r in results[lsq]["data"]}
    n256=f"{row_data.get(256,float('nan')):>6.2f}" if 256 in row_data else "  N/A "
    print(f"{lsq:>6}  {label:<45}  {row_data.get(32,99):>6.2f}  {row_data.get(64,99):>6.2f}  {row_data.get(128,99):>6.2f}  {n256}")

print()
print("=== DIVERGENCE TEST ===")
print("Stable (near-PSD): log|ev0| stays near -11..")
print("Diverging: log|ev0| grows toward 0 (or positive) as N increases")
for lsq,label in TEST_LSQS:
    row_data={r["N"]:r["log10_abs_ev0"] for r in results[lsq]["data"]}
    v32=row_data.get(32,0); v64=row_data.get(64,0); v128=row_data.get(128,0)
    # Diverging if log|ev0| increases by >2 from N=32 to N=128
    delta=v128-v32
    verdict="STABLE" if abs(delta)<2 else "DIVERGING"
    print(f"  lsq={lsq:>4}: Δlog(N=32→128)={delta:>+6.2f}  → {verdict}")

os.makedirs("results",exist_ok=True)
with open("results/N_scaling_v2.json","w") as f:
    json.dump({"results":{str(lsq):v for lsq,v in results.items()}},f,indent=2)
print(f"\nSaved → results/N_scaling_v2.json")

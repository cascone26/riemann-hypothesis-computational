#!/usr/bin/env python3
"""
small_j_hunt.py — Test γ₅, γ₈, γ₉ at smaller j values (lsq < 450).

Prior hunts tested lsq=[450,1500]. This covers lsq=[250,470].
If these zeros have stable resonances at smaller lsq we'd have missed them.

γ₅=32.935: j=29→256, j=30→307, j=31→372  (j=32→450 pseudo-stable, j=33→542 pseudo-stable)
γ₈=43.327: j=40→329, j=41→390, j=42→468  (j=43→508 pseudo-stable)
γ₉=48.005: j=44→318?, j=45→360, j=46→406, j=47→471  (j=48→535 diverge)

Also test lsq in {316,318,319,320} for γ₉ (close to j=44 natural lsq).
Also: small j for γ₅ at j=28 → lsq~210 (borderline on having enough primes).
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900
MAX_LSQ=500
sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sv[i]: sv[i*i::i]=False
all_primes=np.where(sv)[0]

ZEROS={5:32.935061587739189, 8:43.327073280914999, 9:48.005150881167159}
print("Natural j→lsq mapping for untested region:")
for k,g in ZEROS.items():
    for j in range(26,50):
        lsq=round(math.exp(TWO_PI*j/g))
        if 200<=lsq<=480:
            print(f"  γ_{k}(g={g:.3f}) j={j} → lsq={lsq}")
print()

TESTS=[
    (5,  ZEROS[5], 28,  None),  # γ₅ j=28
    (5,  ZEROS[5], 29,  None),  # γ₅ j=29
    (5,  ZEROS[5], 30,  None),  # γ₅ j=30
    (5,  ZEROS[5], 31,  None),  # γ₅ j=31
    (8,  ZEROS[8], 39,  None),  # γ₈ j=39
    (8,  ZEROS[8], 40,  None),  # γ₈ j=40
    (8,  ZEROS[8], 41,  None),  # γ₈ j=41
    (8,  ZEROS[8], 42,  None),  # γ₈ j=42
    (9,  ZEROS[9], 44,  None),  # γ₉ j=44
    (9,  ZEROS[9], 45,  None),  # γ₉ j=45
    (9,  ZEROS[9], 46,  None),  # γ₉ j=46
    (9,  ZEROS[9], 47,  None),  # γ₉ j=47
]
# Fill in lsq
TESTS=[(k,g,j,round(math.exp(TWO_PI*j/g))) for k,g,j,_ in TESTS]

N_VALUES=[32,64,96,128,192]

def compute_ev0(lsq, N):
    if lsq>MAX_LSQ: return float('nan')
    G=G_PER_N*N; L=2*math.log(math.sqrt(lsq)); ns=np.arange(-N,N+1); D=len(ns)
    R=TWO_PI*N/L+50.0; grid=np.linspace(-R,R,G)
    h=grid[1]-grid[0]; w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=np.real(digamma(0.25+0.5j*grid))-LOGPI
    c=TWO_PI*ns[:,None]/L; d=c-grid[None,:]
    sm=np.abs(d)<1e-12
    E=np.where(sm,np.sqrt(L)+0j,(1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch=(E*(w*Om)@E.conj().T).real/TWO_PI
    sh=np.sinh(L/4)**2; ni=ns[:,None]; nj=ns[None,:]
    B=32*L*sh*(L**2-16*math.pi**2*ni*nj)/((L**2+16*math.pi**2*ni**2)*(L**2+16*math.pi**2*nj**2))
    p_in=all_primes[all_primes<=lsq]; ll=[]; ww=[]
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

results={}
for k,g,j,lsq in TESTS:
    if lsq>MAX_LSQ:
        print(f"SKIP γ_{k} j={j} lsq={lsq} (>MAX_LSQ={MAX_LSQ})")
        continue
    print(f"\n=== γ_{k}(γ={g:.3f}) j={j} → lsq={lsq} ===", flush=True)
    row=[]
    for N in N_VALUES:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"log10_abs_ev0":log_e})
        marker=" NEAR-PSD" if log_e<-10 else ""
        print(f"  N={N:>3}: log={log_e:>7.2f}  {dt:.1f}s{marker}", flush=True)
    results[f"k{k}_j{j}_lsq{lsq}"]={"k":k,"j":j,"lsq":lsq,"data":row}

print("\n\n=== SUMMARY ===")
print(f"{'zero':>10}  {'j':>4}  {'lsq':>5}  N32  N64  N128 N192  → verdict")
print("-"*75)
for k,g,j,lsq in TESTS:
    key=f"k{k}_j{j}_lsq{lsq}"
    if key not in results: continue
    rd={r["N"]:r["log10_abs_ev0"] for r in results[key]["data"]}
    delta=rd[192]-rd[32] if 192 in rd and 32 in rd else 99
    maxlog=max(rd.values())
    stable=abs(delta)<2 and maxlog<-10
    pseudo=rd.get(128,-99)<-10 and not stable
    verdict="STABLE ***" if stable else ("PSEUDO" if pseudo else "DIVERGING")
    vals=" ".join(f"{rd.get(N,-99):>5.1f}" for N in [32,64,128,192])
    print(f"γ_{k:>2}(j={j:>2})  {lsq:>5}  {vals}  → {verdict}")

os.makedirs("results",exist_ok=True)
with open("results/small_j_hunt.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/small_j_hunt.json")

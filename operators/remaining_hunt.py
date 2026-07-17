#!/usr/bin/env python3
"""
remaining_hunt.py — Finish the hunt for γ₅, γ₈, γ₉ stable resonances.

Untested j values for these zeros:
- γ₅=32.935: j=34→lsq=653 (j=33/542 diverged at N=192; test j=34 + re-check 542 at N=256)
- γ₈=43.327: j=43→508, j=44→591, j=45→683 (j=46/789 and j=47/912 diverged)
- γ₉=48.005: j=48→535, j=49→609, j=50→697 (j=51/792 and j=52/903 diverged)

Also: γ₅ j=33/lsq=542 extended to N=256 to check if it recovers like lsq=797 did.
Also: lsq=797 extended to N=512,640 to determine ultimate stability.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

# γ₅ j=33 extended, + untested j values
# γ₈ untested j values
# γ₉ untested j values
# lsq=797 extended
TESTS=[
    (5,  32.935061587739189, 33,  542, [192, 256]),         # γ₅ extend
    (5,  32.935061587739189, 34,  653, [32,64,96,128,192]), # γ₅ new
    (8,  43.327073280914999, 43,  508, [32,64,96,128,192]),
    (8,  43.327073280914999, 44,  591, [32,64,96,128,192]),
    (8,  43.327073280914999, 45,  683, [32,64,96,128,192]),
    (9,  48.005150881167159, 48,  535, [32,64,96,128,192]),
    (9,  48.005150881167159, 49,  609, [32,64,96,128,192]),
    (9,  48.005150881167159, 50,  697, [32,64,96,128,192]),
    (0,  0.0,                 0,  797, [448,512,640]),       # lsq=797 mystery extend
]

MAX_LSQ=750
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
for k,gamma,j,lsq,N_vals in TESTS:
    label=f"γ_{k}={gamma:.3f} j={j}" if k>0 else f"lsq={lsq} mystery"
    print(f"\n=== {label} → lsq={lsq} ===")
    row=[]
    all_logs=[]
    for N in N_vals:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        all_logs.append(log_e)
        row.append({"N":N,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
        marker=" NEAR-PSD" if log_e<-10 else ""
        print(f"  N={N:>3}: ev0={e:>12.4e}  log={log_e:>7.2f}  {dt:.1f}s{marker}", flush=True)
    key=f"k{k}_j{j}_lsq{lsq}"
    results[key]={"k":k,"gamma":gamma,"j":j,"lsq":lsq,"data":row}

print("\n\n=== SUMMARY ===")
print(f"{'zero':>10}  {'j':>4}  {'lsq':>5}  logs_at_N  verdict")
print("-"*75)
for k,gamma,j,lsq,N_vals in TESTS:
    key=f"k{k}_j{j}_lsq{lsq}"
    if key not in results: continue
    rd={r["N"]:r["log10_abs_ev0"] for r in results[key]["data"]}
    nmin=min(rd.keys()); nmax=max(rd.keys())
    delta=rd[nmax]-rd[nmin]
    verdict="STABLE" if abs(delta)<2 else "DIVERGING"
    vals=" ".join(f"N{N}:{rd[N]:.1f}" for N in sorted(rd))
    label=f"γ_{k}(j={j})" if k>0 else f"lsq=797"
    star=" ***" if verdict=="STABLE" else ""
    print(f"{label:>12}  {lsq:>5}  {vals}  → {verdict}{star}")

os.makedirs("results",exist_ok=True)
with open("results/remaining_hunt.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/remaining_hunt.json")

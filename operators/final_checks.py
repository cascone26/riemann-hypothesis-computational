#!/usr/bin/env python3
"""
final_checks.py — Targeted high-N checks.

1. lsq=797 at N=512,640 (correct MAX_LSQ=1000 sieve — no prime omission)
2. lsq=967 at N=320,384 (confirm γ₂ stable at high N, matching lsq=787 pattern)
3. γ₈ j=43 (lsq=508) extended to N=256 (pseudo-stable to N=128, does it recover?)
4. γ₅ j=32 (lsq=450) — untested j just below the j=33/lsq=542 pseudo-stable one

Also: endpoint vanishing check — verify that prime lsq contributes S_n=0.
"""
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

# Endpoint vanishing check
print("=== ENDPOINT VANISHING CHECK ===")
for lsq in [787, 797, 967]:
    L=math.log(lsq)
    # Check if lsq is prime
    is_prime = lsq<=MAX_LSQ and sv[lsq]
    if is_prime and lsq<=MAX_LSQ:
        # Compute S_1 contribution from p=lsq
        contrib_S1 = (math.log(lsq)/math.sqrt(lsq)) * math.sin(2*math.pi * 1 * math.log(lsq)/L)
        contrib_Sd0 = (math.log(lsq)/math.sqrt(lsq)) * 2 * math.cos(0) * (1 - math.log(lsq)/L)
        print(f"  lsq={lsq} (prime): S_1 contrib={contrib_S1:.2e}, Sd_0 contrib={contrib_Sd0:.2e}")
print()

TESTS=[
    ("lsq=797 ext",  0, 0,  797, [512, 640]),
    ("lsq=967 ext", 2, 23, 967, [320, 384]),
    ("γ₈ j=43 ext", 8, 43, 508, [192, 256]),
    ("γ₅ j=32 new", 5, 32, 450, [32, 64, 96, 128, 192]),
]

results={}
for label,k,j,lsq,N_vals in TESTS:
    print(f"\n=== {label}: lsq={lsq} ===")
    row=[]
    for N in N_vals:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
        marker=" STABLE" if log_e<-10 else " DIVERGED"
        print(f"  N={N:>3}: ev0={e:>12.4e}  log={log_e:>7.2f}  {dt:.1f}s{marker}", flush=True)
    results[label]={"lsq":lsq,"data":row}

print("\n\n=== FINAL SUMMARY ===")
for label,k,j,lsq,N_vals in TESTS:
    rd={r["N"]:r["log10_abs_ev0"] for r in results[label]["data"]}
    nmin=min(rd.keys()); nmax=max(rd.keys())
    delta=rd[nmax]-rd[nmin]
    maxlog=max(rd.values())
    stable_near_psd=(abs(delta)<2) and (maxlog<-10)
    verdict="STABLE+NEARPSD" if stable_near_psd else ("STABLE_NOT_PSD" if abs(delta)<2 else "DIVERGING")
    vals=" ".join(f"N{N}:{rd[N]:.2f}" for N in sorted(rd))
    print(f"  {label}: {vals} → {verdict}")

os.makedirs("results",exist_ok=True)
with open("results/final_checks.json","w") as f:
    json.dump(results,f,indent=2)
print("\nSaved → results/final_checks.json")

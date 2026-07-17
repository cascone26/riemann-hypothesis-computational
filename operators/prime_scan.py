#!/usr/bin/env python3
"""
prime_scan.py — Test ALL primes in [700,810] for N-scaling stability.

CRITICAL QUESTION: Is stability specific to Riemann-resonant primes (787, 797)
or is every prime lsq in this range stably near-PSD?

If ALL primes are stable: endpoint vanishing makes ALL prime lsq near-PSD;
the N-scaling test doesn't discriminate; our main result needs rethinking.

If ONLY resonant primes are stable: zero resonance is genuinely selective;
endpoint vanishing is a necessary but not sufficient condition.

Known resonant primes in [700,810]:
- lsq=787: γ₁ resonance (j=15)
- lsq=797: mystery (no zero match found)
- lsq=709: γ₁₀ resonance (j=52)

All other primes in range should diverge if zero resonance is the key mechanism.
"""
import math, numpy as np, scipy.linalg, time, json, os
from scipy.special import digamma

TWO_PI=2*math.pi; LOGPI=math.log(math.pi)
G_PER_N=900

MAX_LSQ=820
sv=np.ones(MAX_LSQ+1,dtype=bool); sv[0]=sv[1]=False
for i in range(2,int(MAX_LSQ**0.5)+1):
    if sv[i]: sv[i*i::i]=False
primes_in_range=[p for p in range(700,811) if sv[p]]
all_primes=np.where(sv)[0]

print(f"Primes in [700,810]: {primes_in_range}")
print(f"Total: {len(primes_in_range)} primes")
print()

# Known resonances in range
resonances = {
    709: "γ₁₀ resonance (j=52)",
    787: "γ₁ resonance (j=15)",
    797: "mystery stable",
}
for p,note in resonances.items():
    print(f"  p={p}: {note}")
print()

N_VALUES=[32, 64, 128, 192]

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
for lsq in primes_in_range:
    row=[]
    for N in N_VALUES:
        t0=time.time()
        e=compute_ev0(lsq,N)
        dt=time.time()-t0
        log_e=math.log10(abs(e)) if abs(e)>0 else -99
        row.append({"N":N,"ev0":e,"log10_abs_ev0":log_e,"time_s":dt})
    results[lsq]=row
    rd={r["N"]:r["log10_abs_ev0"] for r in row}
    delta=rd[192]-rd[32]
    stable=abs(delta)<2 and rd[192]<-10
    note=f"  ← {resonances[lsq]}" if lsq in resonances else ""
    verdict="STABLE" if stable else ("PSEUDO" if rd[128]<-10 and rd[192]>-10 else "DIVERGING")
    print(f"p={lsq:>4}:  N32={rd[32]:>6.2f}  N64={rd[64]:>6.2f}  N128={rd[128]:>6.2f}  N192={rd[192]:>6.2f}  → {verdict}{note}", flush=True)

print()
print("=== VERDICT ===")
stable_primes=[]; diverging_primes=[]; pseudo_primes=[]
for lsq in primes_in_range:
    rd={r["N"]:r["log10_abs_ev0"] for r in results[lsq]}
    delta=rd[192]-rd[32]; stable=abs(delta)<2 and rd[192]<-10
    pseudo=rd[128]<-10 and not stable
    if stable: stable_primes.append(lsq)
    elif pseudo: pseudo_primes.append(lsq)
    else: diverging_primes.append(lsq)

print(f"STABLE prime lsq (N-scaling stable, near-PSD at N=192): {stable_primes}")
print(f"PSEUDO prime lsq (near-PSD to N=128, diverges at N=192): {pseudo_primes}")
print(f"DIVERGING prime lsq: {diverging_primes}")
print()
print("KEY FINDING:")
if len(stable_primes)>5:
    print("MOST primes are stable → endpoint vanishing creates broad prime near-PSD; zero resonance NOT selective at this N.")
elif set(stable_primes).issubset({709,787,797}):
    print("ONLY resonant primes (709,787,797) are stable → zero resonance IS selective even among primes!")
else:
    print(f"Intermediate: {len(stable_primes)}/{len(primes_in_range)} primes stable. Need deeper analysis.")

os.makedirs("results",exist_ok=True)
with open("results/prime_scan.json","w") as f:
    json.dump({str(k):v for k,v in results.items()},f,indent=2)
print("Saved → results/prime_scan.json")

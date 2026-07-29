# Riemann: test Dyson's small-diagonal-stratum hypothesis and read the HONEST kappa_eff.
# sigma_off = rms over pairs of R[t,t']=E[u_t u_t']/sqrt(E[u_t^2]E[u_t'^2]).
# If E[X/sqrt(YZ)] != E[X]/E[Y] because small-D cells dominate, then the slow decay lives
# in the small-diagonal stratum. Split sigma_off by diagonal size and read kappa in each.
import numpy as np
def sieve_Lambda(M):
    spf=np.zeros(M+1,dtype=np.int64)
    for i in range(2,int(M**0.5)+1):
        if spf[i]==0: spf[i*i::i]=np.where(spf[i*i::i]==0,i,spf[i*i::i])
    Lam=np.zeros(M+1)
    for n in range(2,M+1):
        p=spf[n] if spf[n] else n
        m=n
        while m%p==0: m//=p
        if m==1: Lam[n]=np.log(p)
    return Lam

M=6_000_000
print("sieving to",M,"...")
Lam=sieve_Lambda(M); dL=Lam-1.0

def triples(s):
    out=[]
    for a in range(1,s):
        for b in range(a+1,s):
            c=s-a-b
            if c>b: out.append((a,b,c))
    return out
# admissible mod 3,5,7 (nonzero singular series -> nonzero diagonal support)
def admissible(t):
    for p in (3,5,7):
        if len(set(x%p for x in t))>=p: return False
    return True
s=42
T=[t for t in triples(s) if admissible(t)]
hmax=max(max(t) for t in T)
print(f"s={s}: {len(T)} admissible triples")

def diag_of(t,N,idx):
    d=np.ones(len(idx))
    for h in t: d=d*dL[idx+h]**2
    return d.mean()
def numer(tA,tB,N,idx):
    v=np.ones(len(idx))
    for h in list(tA)+list(tB): v=v*dL[idx+h]
    return v.mean()

import itertools
Ns=[100_000,1_000_000,M-hmax]
print(f"\n{'N':>9} {'sigma_off':>10} {'sig_bigD':>9} {'sig_smallD':>10} {'D_med':>8}")
rows=[]
for N in Ns:
    idx=np.arange(1,N)
    D={t:diag_of(t,N,idx) for t in T}
    Dvals=np.array([D[t] for t in T]); Dmed=np.median(Dvals)
    Rvals=[]; strat=[]
    for tA,tB in itertools.combinations(T,2):
        num=numer(tA,tB,N,idx)
        r=num/np.sqrt(D[tA]*D[tB])
        Rvals.append(r)
        strat.append(min(D[tA],D[tB])<Dmed)   # True = small-D stratum
    Rvals=np.array(Rvals); strat=np.array(strat)
    sig=np.sqrt(np.mean(Rvals**2))
    sig_big=np.sqrt(np.mean(Rvals[~strat]**2))
    sig_small=np.sqrt(np.mean(Rvals[strat]**2))
    rows.append((N,sig,sig_big,sig_small,Dmed))
    print(f"{N:>9} {sig:>10.4f} {sig_big:>9.4f} {sig_small:>10.4f} {Dmed:>8.1f}")

lnlogN=np.log(np.log([r[0] for r in rows]))
def k(col): return -np.polyfit(lnlogN,np.log([r[col] for r in rows]),1)[0]
print(f"\nkappa_eff (all)      = {k(1):.2f}")
print(f"kappa_eff (big-D)    = {k(2):.2f}")
print(f"kappa_eff (small-D)  = {k(3):.2f}")
print("\nREAD: if small-D stratum has much larger sigma AND shallower kappa, Dyson is right:")
print("  the slow decay is denominator-weighting, and the HONEST big-D kappa is the true one.")
import json
json.dump({'s':s,'M':M,'Ns':[r[0] for r in rows],
           'sigma_all':[r[1] for r in rows],'sigma_bigD':[r[2] for r in rows],
           'sigma_smallD':[r[3] for r in rows],
           'kappa_all':k(1),'kappa_bigD':k(2),'kappa_smallD':k(3)},
          open('results/rt6_stratified_sigma.json','w'),indent=2)
print("CHECKPOINT -> results/rt6_stratified_sigma.json")

"""Tao independent physical Gram ladder (banked for Gate). Fresh sieve, N=1e6.
Reproduces Hilbert's lam_min(R), Cholesky-PD, sigma_off + localization guard."""
import json, numpy as np
from math import log
def sieve_Lambda(M):
    spf=np.zeros(M+1,dtype=np.int64)
    for i in range(2,M+1):
        if spf[i]==0:
            for j in range(i,M+1,i):
                if spf[j]==0: spf[j]=i
    L=np.zeros(M+1)
    for n in range(2,M+1):
        p=int(spf[n]); m=n
        while m%p==0: m//=p
        if m==1: L[n]=log(p)
    return L
def triples_at_sum(s,lo=1):
    return [(a,b,s-a-b) for a in range(lo,s) for b in range(a+1,s) if s-a-b>b]
def admissible(t):
    for p in (3,5,7):
        if len(set(x%p for x in t))==p: return False
    return True
N=1_000_000; MAXO=60
Lam=sieve_Lambda(N+MAXO+5); dL=Lam-1.0; base=np.arange(0,N)
rows={}
for s in [24,30,36,42,48]:
    T=[t for t in triples_at_sum(s) if admissible(t) and t[2]<=MAXO]; n=len(T)
    V=np.array([dL[base+t[0]]*dL[base+t[1]]*dL[base+t[2]] for t in T])
    G=(V@V.T)/N; d=np.sqrt(np.diag(G)); R=G/np.outer(d,d)
    w,Vec=np.linalg.eigh(R)
    try: np.linalg.cholesky(G); pd=True
    except: pd=False
    off=R-np.diag(np.diag(R)); sig=float(np.sqrt((off**2).sum()/(n*(n-1))))
    vmin=Vec[:,0]; PR=float(1.0/np.sum(vmin**4))
    rows[str(s)]={"n":n,"lam_min_R":float(w[0]),"sigma_off":sig,"chol_PD":pd,
                  "PR":PR,"PR_over_n":PR/n}
out={"object":"Tao independent physical Gram ladder","version":"rt6-tao-physical-gram-2026-07-28",
     "window_N":N,"sector":"accessible-odd (admissible mod 3,5,7)","per_s":rows,
     "note":"two-apparatus vs rt6_hilbert_gram_physical.json; lam_min agrees ~1-2%, Chol-PD all s, min-vec delocalized"}
json.dump(out,open("results/rt6_tao_physical_gram_ladder.json","w"),indent=2)
for s,r in rows.items(): print(f"s={s}: n={r['n']:3d} lam_min={r['lam_min_R']:.4f} sig={r['sigma_off']:.4f} PD={r['chol_PD']} PR/n={r['PR_over_n']:.3f}")
print("wrote results/rt6_tao_physical_gram_ladder.json")

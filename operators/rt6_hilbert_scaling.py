import numpy as np
def build_dL(N):
    LIM=N+70
    spf=np.zeros(LIM,dtype=np.int64)
    for i in range(2,LIM):
        if spf[i]==0: spf[i:LIM:i]=np.where(spf[i:LIM:i]==0,i,spf[i:LIM:i])
    Lam=np.zeros(LIM)
    for n in range(2,LIM):
        p=int(spf[n]); m=n
        while m%p==0: m//=p
        if m==1: Lam[n]=np.log(p)
    return Lam-1.0

def admissible(t):
    for p in (3,5,7):
        if len({x%p for x in t})==p: return False
    return True
def triples(s):
    return [(a,b,c) for a in range(1,s) for b in range(a+1,s)
            for c in (s-a-b,) if c>b and admissible((a,b,c))]
def stats(dL,s,base,span):
    T=triples(s); n=np.arange(base,base+span)
    U=np.stack([dL[n+t[0]]*dL[n+t[1]]*dL[n+t[2]] for t in T])
    G=(U@U.T)/span; d=np.sqrt(np.diag(G)); R=G/np.outer(d,d)
    lam=np.linalg.eigvalsh(R).min()
    off=R-np.eye(len(T)); so=np.sqrt((off**2).sum()/(len(T)*(len(T)-1)))
    return len(T),lam,so

N=4_000_000
print("building dL, N=4e6 ...",flush=True); dL=build_dL(N)

print("\n=== AXIS 1: ||O_norm||=1-lam_min(R) vs s at FIXED N=4e6 (pin the exponent) ===")
rows=[]
for s in (24,30,36,42,48,54,60,66,72):
    n,lam,so=stats(dL,s,2,N-80); Onorm=1-lam
    rows.append((s,n,Onorm,so)); print(f"  s={s:2d} n={n:3d}  ||O||={Onorm:.4f}  sigma_off={so:.4f}")
import numpy as np
ns=np.array([r[1] for r in rows]); On=np.array([r[2] for r in rows])
sl_all=np.polyfit(np.log(ns),np.log(On),1)[0]
sl_drop48=np.polyfit(np.log([r[1] for r in rows if r[0]!=48]),np.log([r[2] for r in rows if r[0]!=48]),1)[0]
print(f"  log-log slope ||O|| vs n: ALL={sl_all:.3f}  (drop s=48)={sl_drop48:.3f}   [generic sqrt = 0.5]")

print("\n=== AXIS 2: sigma_off vs N at FIXED s=42 (pin the decay alpha) ===")
srow=[]
for span in (250000,500000,1000000,2000000,3900000):
    n,lam,so=stats(dL,42,2,span); srow.append((span,so)); print(f"  N={span:>8d}  sigma_off={so:.5f}  lam_min={lam:.4f}")
Ns=np.array([r[0] for r in srow]); ss=np.array([r[1] for r in srow])
alpha=-np.polyfit(np.log(Ns),np.log(ss),1)[0]
print(f"  sigma_off ~ N^(-alpha), alpha={alpha:.3f}")
print(f"\n  COUPLING TEST: sigma_off*sqrt(n) at fixed N=4e6 vs s: ",
      [f'{r[3]*np.sqrt(r[1]):.3f}' for r in rows])

# Riemann -> Hilbert: the WORST-conditioned 5x5 admissible (mod-105) principal cell,
# so the Arb minor certifies the tightest leading minor, not a random one.
# Relative off-diagonal structure only (scale-free) = my kernel's validated authority.
from fractions import Fraction as Fr
from math import factorial
from itertools import combinations
import numpy as np
ODD=[3,5,7,11,13]
def wp(pos,p):
    d=set(pos); k=len(d)
    if k==0: return Fr(1)
    nu=len(set(x%p for x in d))
    return Fr(p-nu,p)/Fr(p-1,p)**k
def W(cfg):
    v=Fr(1)
    for p in ODD: v*=wp(cfg,p)
    return v
def sp(c):
    c=list(c)
    if len(c)==1: yield [c]; return
    f=c[0]
    for s in sp(c[1:]):
        for i,x in enumerate(s): yield s[:i]+[[f]+x]+s[i+1:]
        yield [[f]]+s
def bmu(n): return ((-1)**(n-1))*factorial(n-1)
PARTS=[]
for pi in sp(range(6)):
    mu=1
    for B in pi: mu*=bmu(len(B))
    PARTS.append((mu,[tuple(B) for B in pi]))
def O(X):
    t=Fr(0)
    for mu,pi in PARTS:
        pr=Fr(mu)
        for B in pi: pr*=W([X[i] for i in B])
        t+=pr
    return t
def triples(s):
    out=[]
    for a in range(1,s):
        for b in range(a+1,s):
            c=s-a-b
            if c>b: out.append((a,b,c))
    return out

s=48
T=[t for t in triples(s) if float(W(list(t)))>0]     # admissible-odd sector
n=len(T)
# scale-free relative R: normalize off by sqrt(Wself Wself); diagonal set to per-row max scale
Ws=[float(W(list(t))) for t in T]
Rmat=np.eye(n)
for i in range(n):
    for j in range(i+1,n):
        o=float(O(list(T[i])+list(T[j])))
        r=o/np.sqrt(Ws[i]*Ws[j])
        Rmat[i,j]=r; Rmat[j,i]=r
# rank cells by |offdiag| row-strength; take the 12 strongest, exhaustively find worst 5x5
strength=np.sum(np.abs(Rmat-np.eye(n)),axis=1)
cand=list(np.argsort(strength)[::-1][:12])
best=None
for combo in combinations(cand,5):
    sub=Rmat[np.ix_(combo,combo)]
    # relative conditioning: ratio of off-Frobenius to diag -> proxy for how hard PD is
    off=np.sqrt(np.sum((sub-np.diag(np.diag(sub)))**2))
    lam=np.linalg.eigvalsh(sub).min()   # relative (scale-free) min eig of the 5x5 structure
    score=off                            # largest off-diagonal mass = hardest cell
    if best is None or score>best[0]:
        best=(score,combo,lam,off)
_,combo,lam,off=best
print(f"s={s}, admissible-odd sector n={n}")
print("WORST-CONDITIONED 5x5 admissible cell (largest off-diagonal mass) -> Hilbert's minor target:")
for idx in combo:
    print(f"   t={T[idx]}   (mod3 res {sorted(set(x%3 for x in T[idx]))}, "
          f"mod5 {sorted(set(x%5 for x in T[idx]))}, mod7 {sorted(set(x%7 for x in T[idx]))})")
print(f"   relative off-diag Frobenius mass = {off:.3f}  (scale-free; the tight leading minor)")
print("   NOTE: values are relative structure only; Hilbert reads physical magnitudes from his Gram G.")

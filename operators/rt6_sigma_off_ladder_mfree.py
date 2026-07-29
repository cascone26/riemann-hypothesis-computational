# Dyson's all-orders decider (Riemann's lane): does sigma_off(s)*sqrt(n) stay < 1,
# or does it grow like the generic bulk-edge (pure dimension effect)?
# R[t,t'] = O[t,t'] / sqrt(D_t D_t'),  D_t = E[u_t^2] = m2^3 * Wself(t)  (SQUARED diag, >=0).
# O = connected kernel off-diagonal (3-apparatus confirmed).  m2 = E[dLambda^2] ~ log N.
from fractions import Fraction as Fr
import numpy as np, json
ODD=[3,5,7,11,13]
M2 = 13.20   # Dyson's E[dLambda^2] at N=4e6 (log-enhanced; common scale)
def wp(pos,p):
    d=set(pos); k=len(d)
    if k==0: return Fr(1)
    nu=len(set(x%p for x in d))
    return Fr(p-nu,p)/Fr(p-1,p)**k
def W(cfg):
    v=Fr(1)
    for p in ODD: v*=wp(cfg,p)
    return v
from math import factorial
def set_parts(c):
    c=list(c)
    if len(c)==1: yield [c]; return
    f=c[0]
    for sm in set_parts(c[1:]):
        for i,s in enumerate(sm): yield sm[:i]+[[f]+s]+sm[i+1:]
        yield [[f]]+sm
def bmu(n): return ((-1)**(n-1))*factorial(n-1)
PARTS6=[]
for pi in set_parts(range(6)):
    mu=1
    for B in pi: mu*=bmu(len(B))
    PARTS6.append((mu,[tuple(B) for B in pi]))
def O_conn(X):
    tot=Fr(0)
    for mu,pi in PARTS6:
        pr=Fr(mu)
        for B in pi: pr*=W([X[i] for i in B])
        tot+=pr
    return tot
def triples(s):
    out=[]
    for a in range(1,s):
        for b in range(a+1,s):
            c=s-a-b
            if c>b: out.append((a,b,c))
    return out

rows=[]
for s in (30,36,42,48,54):
    T0=triples(s)
    T=[t for t in T0 if float(W(list(t)))>0]   # accessible odd sector: admissible triples only
    n=len(T)
    D=np.array([float(W(list(t))) for t in T])   # m2-FREE: field scale cancels in R (consistent units)
    sq=np.sqrt(D)
    R=np.zeros((n,n))
    for i in range(n):
        R[i,i]=1.0
        for j in range(i+1,n):
            o=float(O_conn(list(T[i])+list(T[j])))
            r=o/(sq[i]*sq[j])
            R[i,j]=r; R[j,i]=r
    ev=np.linalg.eigvalsh(R)
    off=R[~np.eye(n,dtype=bool)]
    sigma=float(np.sqrt(np.mean(off**2)))      # rms off-diagonal entry
    prod=sigma*np.sqrt(n)
    edge=1.0-2.0*prod                          # Wigner bulk-edge prediction for lambda_min
    rows.append((s,n,float(ev.min()),sigma,prod,edge,float(D.min()),float(D.max())))
    print(f"s={s:2d} n={n:3d}  lmin(R)={ev.min():+.4f}  sigma_off={sigma:.4f}  "
          f"sigma*sqrt(n)={prod:.4f}  bulk-edge(1-2p)={edge:+.4f}")

print("\n=== ALL-ORDERS DECIDER: is sigma_off*sqrt(n) bounded < 0.5 (=> lmin>0), or -> 0.5+ ? ===")
ps=[r[4] for r in rows]
print("  sigma*sqrt(n) ladder:", "  ".join(f"{p:.4f}" for p in ps))
print("  lmin(R) vs bulk-edge (1-2p): does arithmetic beat the generic null?")
for r in rows:
    gap=r[2]-r[5]
    print(f"    s={r[0]}: lmin={r[2]:+.4f}  edge={r[5]:+.4f}  arith_residue={gap:+.4f} "
          f"({'ABOVE null (arith helps)' if gap>0.02 else 'tracks generic null' if abs(gap)<=0.02 else 'BELOW null'})")
json.dump([{'s':r[0],'n':r[1],'lmin':r[2],'sigma_off':r[3],'sigma_sqrtn':r[4],
            'bulk_edge':r[5],'D_min':r[6],'D_max':r[7]} for r in rows],
          open('results/rt6_sigma_off_ladder.json','w'), indent=2)
print("\nCHECKPOINT -> results/rt6_sigma_off_ladder.json  (m2=13.20)")

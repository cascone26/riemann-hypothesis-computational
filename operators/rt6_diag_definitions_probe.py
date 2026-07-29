# Completeness form as the COVARIANCE matrix of T(t)=a_w1 a_w2 a_w3:
#   M[t,t'] = W(union config) - W(t)*W(t')   (connected/covariance, PSD by construction)
# W(cfg) = Euler product of collapse-corrected local density over odd primes.
# Diagonal M[t,t] = W(t) - W(t)^2 = W(t)(1-W(t)) > 0 for 0<W(t)<1. No rho, no kappa4.
from fractions import Fraction as Fr
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
def triples(s):
    out=[]
    for a in range(1,s):
        for b in range(a+1,s):
            c=s-a-b
            if c>b: out.append((a,b,c))
    return out
for s in (30,36,42):
    T=triples(s); n=len(T)
    Wt=[W(t) for t in T]
    Mrat=[[W(list(T[i])+list(T[j]))-Wt[i]*Wt[j] for j in range(n)] for i in range(n)]
    M=np.array([[float(x) for x in row] for row in Mrat])
    M=(M+M.T)/2
    ev=np.linalg.eigvalsh(M)
    diag=np.diag(M)
    off=M-np.diag(diag)
    offnorm=float(np.max(np.abs(np.linalg.eigvalsh(off))))
    dmin=float(diag.min()); dmax=float(diag.max())
    ratio=offnorm/dmin if dmin>0 else float('inf')
    # participation ratio of min eigenvector (Tao's alignment test, on M)
    w,vv=np.linalg.eigh(M); vmin=vv[:,0]
    PR=(np.sum(vmin**2)**2)/np.sum(vmin**4)
    print(f"s={s}: n={n}  lambda_min(M)={ev.min():+.6e}  lambda_max={ev.max():.4e}")
    print(f"       diag in [{dmin:.4e},{dmax:.4e}] (all>0: {bool((diag>0).all())})  "
          f"||off||/diag_min={ratio:.4f} {'<0.32' if ratio<0.32 else '>=0.32'}")
    print(f"       PR(min-eigvec)={PR:.2f}  PR/n={PR/n:.3f}  (Tao: const PR/n = delocalized/generic)")

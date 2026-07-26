# Obligation 3. Derivation: diagonal dterms give D_kk(t) = k^2 b^2 (2a-t)cos(kbt) - kb sin(kbt).
# IBP the cos term once (boundary vanishes: sin(0)=0, sin(2a*kb)=sin(k*pi)=0, and (2a-t) kills t=2a):
#   Q[k,k] = -(kb/a) * INT_0^{2a} H'(t) * (2a-t) * sin(kbt) dt        [exact, same family as v2]
# Then H'(t) = (1/2)log t + phi(t) with phi in BV  ==>  Q[k,k] ~ log(kb) + c0 + c1*cos(kb*log2 + p0) + O(1/k)
import numpy as np
from scipy.special import digamma, zeta as hurwitz
import math

A=0.549; b=math.pi/(2*A); log2=math.log(2); twoa=2*A; wp=log2/math.sqrt(2)
Om0=digamma(0.25)-math.log(math.pi); NL=4000
ng=800001
t=np.linspace(1e-12,twoa,ng)
Hp=-2*np.exp(0.5*t)+2*np.exp(-0.5*t)-Om0/2+wp*(t>log2)
la=2*np.arange(NL)+0.5; lc=0.25/(np.arange(NL)+0.25)**2
for i0 in range(0,NL,100):
    sl=slice(i0,min(i0+100,NL))
    Hp+= -(lc[sl]*la[sl]*np.exp(-np.multiply.outer(t,la[sl]))).sum(axis=1)

# sign of the singularity: H'(t) - (1/2)log t should be bounded near 0
for tt in (1e-6,1e-4,1e-2):
    i=np.searchsorted(t,tt)
    print(f"  H'({tt:g})={Hp[i]:+.4f}   H'-0.5*log(t)={Hp[i]-0.5*math.log(t[i]):+.4f}")

Q=np.load('/Users/scones/tmp/Q51.npy'); ODD=list(range(1,52,2))
print("\n   k     Q[k,k] direct    identity RHS     rel.err")
for k in (15,19,27,29,39,51):
    q=Q[ODD.index(k),ODD.index(k)]
    r=-(k*b/A)*np.trapezoid(Hp*(twoa-t)*np.sin(k*b*t),t)
    print(f"  {k:2d}   {q:.8f}   {r:.8f}   {abs(r-q)/q:.2e}")

# asymptotic fit: Q[k,k] = log(kb) + c0 + c1*cos(kb*log2) + c2*sin(kb*log2), odd k in [15,51]
ks=np.array([k for k in ODD if k>=15],float)
d=np.array([Q[ODD.index(int(k)),ODD.index(int(k))] for k in ks])
th=ks*b*log2
M=np.column_stack([np.ones_like(ks),np.cos(th),np.sin(th)])
resid=d-np.log(ks*b)
coef,*_=np.linalg.lstsq(M,resid,rcond=None)
fit=np.log(ks*b)+M@coef
amp=math.hypot(coef[1],coef[2])
print(f"\nfit: Q[k,k] ~ log(kb) + ({coef[0]:+.4f}) + {amp:.4f}*cos(kb*log2+phase)")
print(f"max fit residual over k in [15,51]: {np.max(np.abs(fit-d)):.4f}")
print(f"implied floor for k>=25: log(25b)+c0-amp = {math.log(25*b)+coef[0]-amp:.4f}  (measured min d0 = {d[[k>=25 for k in ks]].min():.4f})")

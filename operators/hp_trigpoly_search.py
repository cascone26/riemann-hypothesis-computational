import numpy as np, itertools
from numpy.polynomial import chebyshev as C
def coeffs(factors):
    # factors: list of (b, power): (b+x)^p in x=cos t ; return cos-coefficients via chebyshev
    p=np.array([1.0])
    for b,e in factors:
        for _ in range(e): p=np.polynomial.polynomial.polymul(p,[b,1.0])
    return C.poly2cheb(p)
def score(c):
    if len(c)<2 or c[1]<=c[0] or (c<-1e-12).any(): return 0,None
    eps=np.linspace(1e-3,c[1]-c[0]-1e-3,4000); f=eps*(c[1]-c[0]-eps)/(c[0]+eps); i=f.argmax()
    return f[i]/c[1:].sum(), eps[i]
base=coeffs([(1,1),(0.5,2)]); print('5841 shape',base/base[-1]*1, score(base))
best=(0,)
for b in np.linspace(0,1.2,61):
  for d in np.linspace(0,1.2,61):
    for form in ('sq2','lin_sq2'):
      fac=[(b,2),(d,2)] if form=='sq2' else [(1,1),(b,2),(d,2)]
      c=coeffs(fac); s,e=score(c)
      if s>best[0]: best=(s,form,b,d,e,c/c.min() if c.min()>0 else c)
print(best)
for b in np.linspace(0,1.5,151):
    c=coeffs([(1,1),(b,2)]); s,e=score(c)
    if s>best[0]*0.0 and abs(b-0.5)<1e-9: print('b=.5',s)

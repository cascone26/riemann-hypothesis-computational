import numpy as np
from scipy.special import digamma, zeta as hurwitz
import math

A=0.549; b=math.pi/(2*A); log2=math.log(2); twoa=2*A; wp=log2/math.sqrt(2)
Om0=digamma(0.25)-math.log(math.pi); Phi0=hurwitz(2,0.25); NL=4000

ng=400001
t=np.linspace(1e-12,twoa,ng)
Hp=-2*np.exp(0.5*t)+2*np.exp(-0.5*t)-Om0/2+wp*(t>log2)
la=2*np.arange(NL)+0.5; lc=0.25/(np.arange(NL)+0.25)**2
for i0 in range(0,NL,200):          # chunked Lerch derivative sum
    sl=slice(i0,min(i0+200,NL))
    Hp+= -(lc[sl]*la[sl]*np.exp(-np.multiply.outer(t,la[sl]))).sum(axis=1)

def rhs(k,l):
    f=Hp*np.sin((k+l)*b*t/2)*np.sin((k-l)*b*t/2)
    return (4*k*l/(k**2-l**2))*np.trapezoid(f,t)/A

Q=np.load('/Users/scones/tmp/Q51.npy'); ODD=list(range(1,52,2))
print("   (k,l)      Q[k,l] direct     identity RHS      rel.err")
for (k,l) in [(15,13),(15,17),(27,29),(39,37),(39,41),(39,7),(51,47),(51,7),(15,49),(3,51)]:
    q=Q[ODD.index(k),ODD.index(l)]
    r=rhs(k,l)
    print(f"  ({k:2d},{l:2d})   {q: .8e}   {r: .8e}   {abs(r-q)/abs(q):.2e}")

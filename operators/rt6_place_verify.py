#!/usr/bin/env python3
"""§423b VERIFY §423's (n−1,1) signature is real, not a regularization artifact. Vary: n (#primes), the
regularization width W (exp(−(r/W)²)), and integration cutoff. If signature stays (n−1,1) [exactly ONE negative]
across all, the single-anchor Hodge structure in place-space is REAL. If the negative count changes with W/cutoff,
it's an artifact (like §412/414)."""
import numpy as np
from scipy.special import digamma
from scipy import integrate
def Om(r): return np.real(digamma(0.25+0.5j*r))-np.log(np.pi)
def build(primes,W,cut):
    L=np.log(primes); n=len(primes); M=np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            f=lambda r: Om(r)*np.cos(L[i]*r)*np.cos(L[j]*r)*np.exp(-(r/W)**2)
            v,_=integrate.quad(f,0,cut,limit=500); M[i,j]=M[j,i]=2*v
    return M
def sig(M,tol=1e-8):
    e=np.linalg.eigvalsh(M); m=max(1,abs(e).max())
    return int((e>tol*m).sum()),int((abs(e)<=tol*m).sum()),int((e<-tol*m).sum())
ALLP=[2,3,5,7,11,13,17,19,23,29,31,37]
print("§423b robustness of the (n−1,1) place-form signature.")
print(f"{'n':>3} {'W':>5} {'cut':>5} {'signature(+,0,−)':>16}  most-neg-eig")
for n in (6,8,10,12):
    for W in (25.0,40.0,60.0):
        for cut in (100,150):
            M=build(np.array(ALLP[:n]),W,cut)
            e=np.linalg.eigvalsh(M)
            print(f"{n:>3} {W:>5.0f} {cut:>5} {str(sig(M)):>16}  {e[0]:+.4f}",flush=True)
print("\nreading: '# negative' column CONSTANT at 1 across n,W,cut => the single-anchor Hodge structure is REAL")
print("(place-space carries signature (n−1,1) intrinsically). Varies => regularization artifact.")

"""Nyman-Beurling-Baez-Duarte criterion. RH <=> d_N = dist(1, span{rho_{1/n}, n<=N}) -> 0 in
L^2(0,1), rho_{1/n}(x)={1/(nx)}. Compute d_N^2 = 1 - b^T G^{-1} b; check d_N^2 ~ C/log N (RH);
examine optimal coefficients c = G^{-1} b vs Mobius mu(n) (the BD conjecture / sign-theorem link)."""
import numpy as np
from scipy import integrate
from sympy import mobius
def frac(y): return y - np.floor(y)
# high-accuracy inner products via adaptive quad, handling the sawtooth (many jumps near 0)
def rho(n, x): return frac(1.0/(n*x))
def inner(j, k):
    # int_0^1 {1/(jx)}{1/(kx)} dx ; substitute to concentrate accuracy; use many panels
    f = lambda x: rho(j,x)*rho(k,x)
    val,_ = integrate.quad(f, 1e-6, 1, limit=400, epsabs=1e-10, epsrel=1e-9)
    return val
def inner1(n):
    f = lambda x: rho(n,x)
    val,_ = integrate.quad(f, 1e-6, 1, limit=400, epsabs=1e-10, epsrel=1e-9)
    return val
for N in (5, 10, 20, 40):
    G = np.array([[inner(j,k) for k in range(1,N+1)] for j in range(1,N+1)])
    b = np.array([inner1(n) for n in range(1,N+1)])
    c = np.linalg.solve(G + 1e-12*np.eye(N), b)
    dN2 = 1.0 - b @ c   # ||1||^2 = 1
    # BD optimal-coefficient conjecture: c_n ~ mu(n)(1 - log n/log N)
    mus = np.array([int(mobius(n)) for n in range(1,N+1)])
    bd_pred = mus*(1 - np.log(np.arange(1,N+1))/np.log(N))
    corr = np.corrcoef(c, bd_pred)[0,1]
    print(f"N={N:>3}: d_N^2={dN2:.5f}  d_N^2*log N={dN2*np.log(N):.4f}  corr(c, mu(n)(1-logn/logN))={corr:+.3f}", flush=True)
print("=> d_N^2 ~ C/log N (RH-consistent); c aligns with Mobius (BD conj + sign-theorem link)?")
print("DONE")

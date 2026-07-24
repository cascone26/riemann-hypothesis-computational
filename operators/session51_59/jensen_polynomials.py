# JENSEN POLYNOMIALS (Griffin-Ono-Rolen-Zagier 2019) for Xi, correct sequence a_k = b_k/(2k)! (Taylor coeffs).
# J^{d,n}(X)=sum_j C(d,j) a_{n+j} X^j ; hyperbolic (all real roots) <=> LP at (d,n). RH <=> all hyperbolic.
# Rescale X=Y/t1 (t1=a_{n+1}/a_n) so coefficients are O(1) -> stable numpy roots.
import numpy as np
from math import comb, factorial
def Phi(u,nmax=60):
    u=np.atleast_1d(np.asarray(u,float)); n=np.arange(1,nmax+1)[:,None]
    T=(2*np.pi**2*n**4*np.exp(9*u)[None,:]-3*np.pi*n**2*np.exp(5*u)[None,:])*np.exp(-np.pi*n**2*np.exp(4*u)[None,:])
    return T.sum(0)
U=2.3; M=500000; u=np.linspace(1e-8,U,M); du=u[1]-u[0]; phi=Phi(u)
K=26; b=[np.sum(u**(2*k)*phi)*du for k in range(K+1)]
a=[b[k]/factorial(2*k) for k in range(K+1)]
def hyp(d,n):
    t1=a[n+1]/a[n]
    # scaled coeffs c_j = C(d,j) * (a_{n+j}/a_n) / t1^j  (poly in Y), all O(1)
    cj=[comb(d,j)*(a[n+j]/a[n])/t1**j for j in range(d+1)]
    r=np.roots(cj[::-1])
    maxim=np.max(np.abs(r.imag)); maxre=np.max(np.abs(r.real))+1e-30
    return maxim/maxre < 1e-5
print("Jensen hyperbolicity J^{d,n} (a_k=b_k/(2k)! Taylor coeffs). H=hyperbolic(LP holds), .=not:")
for d in [2,3,4,5,6]:
    res=[]; first=None
    for n in range(0,K-d):
        h=hyp(d,n); res.append('H' if h else '.')
        if h and first is None: first=n
    print(f"  d={d}: {''.join(res)}   first-H at n={first}")
print("\nGORZ 2019 THEOREM: for each d, J^{d,n} hyperbolic for all n>=N(d) (all but finitely many). RH = zero exceptions.")
print("We should SEE: d=2 all-H (Turan, proven all n); higher d = finite low-n '.' then all H.")

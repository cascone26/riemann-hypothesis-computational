"""FAST quadratic Dirichlet L-function Hardy Z-function via the approximate functional equation.
For real primitive χ_d (root number +1), the completed Λ(s)=Λ(1-s). Hardy Z:
  Z(t) = 2 Σ_{n=1}^{X} χ(n)/√n · cos(θ(t) - t log n)   (main AFE, X ~ sqrt(D t /(2π))),
where θ(t) = Im log[ (D/π)^{(1/4+it/2)} Γ((1/2+it+a)/2) ] is the phase, a=0 even /1 odd.
This is O(X) floats per eval, ~100x faster than Hurwitz. Good enough to locate low zeros."""
import numpy as np
from sympy import jacobi_symbol
from math import lgamma
import cmath

def kronecker(a,n):
    if n==0: return 1 if abs(a)==1 else 0
    res=1
    while n%2==0:
        n//=2
        if a%2==0: return 0
        if a%8 in (3,5): res=-res
    a%=n
    if n==1: return res
    return res*int(jacobi_symbol(a,n))

def sqfree(n):
    n=abs(n); k=2
    while k*k<=n:
        if n%(k*k)==0: return False
        k+=1
    return n!=0
def is_fundamental(d):
    if d==1: return False
    if d%4==1: return sqfree(d)
    if d%4==0:
        m=d//4; return (m%4 in (2,3)) and sqfree(m)
    return False

def chi_arr(d, X):
    D=abs(d)
    return np.array([kronecker(d, n% D) if n%D!=0 else 0 for n in range(X+1)], dtype=float)

def _logGamma(z):  # complex lgamma via cmath (Lanczos through cmath.lgamma unavailable) -> use mpmath-free
    # Stirling is fine for our |z|; but use a series-safe: use numpy? no complex lgamma. Use reflection+Stirling.
    # good enough for phase: Stirling for Re>0
    return (z-0.5)*cmath.log(z) - z + 0.5*cmath.log(2*cmath.pi) + 1/(12*z) - 1/(360*z**3)

def theta_phase(t, D, a):
    s=0.5+1j*t
    val=( (0.25+1j*t/2)*cmath.log(D/cmath.pi) ) + _logGamma((s+a)/2)
    return val.imag

def Zfun_fast(t, d, chi, even):
    D=abs(d); a=0 if even else 1
    X=len(chi)-1
    n=np.arange(1,X+1)
    th=theta_phase(t,D,a)
    terms= chi[1:]/np.sqrt(n)*np.cos(th - t*np.log(n))
    return 2*terms.sum()

def make_L(d, tmax=40.0):
    D=abs(d)
    even = (kronecker(d, D-1)==1)
    X=max(30, int(np.sqrt(D*tmax/(2*np.pi)))+10)
    chi=chi_arr(d,X)
    return chi, even, X

def low_zeros_fast(d, m=6, tmax=40.0, step=0.1):
    chi,even,X=make_L(d,tmax)
    zs=[]; pt=0.05; pv=Zfun_fast(pt,d,chi,even); t=step
    while t<tmax and len(zs)<m:
        cv=Zfun_fast(t,d,chi,even)
        if pv*cv<0:
            lo,hi=pt,t; plo=pv
            for _ in range(40):
                mm=(lo+hi)/2; vm=Zfun_fast(mm,d,chi,even)
                if plo*vm<0: hi=mm
                else: lo,plo=mm,vm
            zs.append((lo+hi)/2)
        pt,pv,t=t,cv,t+step
    return zs

if __name__=="__main__":
    # validate vs known: d=-3 ->8.039, d=-4->6.021, d=5->6.648
    for d in [-3,-4,5,8,-7,13]:
        z=low_zeros_fast(d,m=3)
        print(f"  d={d:4d} even={kronecker(d,abs(d)-1)==1}  low zeros={[round(x,3) for x in z]}")

import numpy as np
from scipy.linalg import eigh
from scipy.integrate import quad
# Weighted-Chebyshev spectral solver for lambda_min of Q(v)=INT INT |x-y|log|x-y| v'(x)v'(y),
# w=v'=sum a_n T_n(x)/sqrt(1-x^2), n=1..K.  K_mn via Gauss-Chebyshev (exact for 1/sqrt weight),
# M_mn closed: V_m(cos phi)=-sin(m phi)/m, M_mn=(1/mn) INT_0^pi sin(m phi)sin(n phi)sin phi dphi.
def Mmn(m,n):
    f=lambda p: np.sin(m*p)*np.sin(n*p)*np.sin(p)
    val,_=quad(f,0,np.pi); return val/(m*n)
def solve(K,N):
    # Gauss-Chebyshev nodes/weights for INT g/sqrt(1-x^2) dx ~ (pi/N) sum g(x_i)
    i=np.arange(1,N+1); xi=np.cos((i-0.5)*np.pi/N); wq=np.pi/N
    # T_m(x_i) for m=1..K
    th=np.arccos(xi); Tm=np.cos(np.outer(np.arange(1,K+1),th))   # (K,N)
    # kernel matrix on nodes
    D=np.subtract.outer(xi,xi); aD=np.abs(D); Kk=np.where(aD>1e-300, aD*np.log(aD), 0.0)  # (N,N)
    # K_mn = wq^2 * sum_ij Kk[i,j] Tm[m,i] Tm[n,j] = wq^2 * Tm @ Kk @ Tm.T
    Kmat=(wq*wq)*(Tm@Kk@Tm.T); Kmat=0.5*(Kmat+Kmat.T)
    M=np.array([[Mmn(m,n) for n in range(1,K+1)] for m in range(1,K+1)]); M=0.5*(M+M.T)
    ev=eigh(Kmat,M,eigvals_only=True)
    return ev[0]
print("Weighted-Chebyshev spectral lambda_min (Gauss-Cheb K, closed-form M):")
for K in [8,12,16,24,32,48]:
    N=max(400,8*K)
    lm=solve(K,N); print(f"  K={K:3d} N={N}: lambda={lm:.10f}  mu1={1+0.5*lm:.10f}",flush=True)
print("Compare sine-basis extrap: lambda~-1.5334, mu1~0.2333. Exponential convergence => many more digits.")

# NEW METHOD (out of the positivity basin): NYMAN-BEURLING / BAEZ-DUARTE approximation criterion.
# In L^2(0,1): rho_theta(x) = {theta/x} - theta*{1/x}.  RH <=> dist(1, span{rho_{1/k}: k=1..N}) -> 0.
# Compute d_N^2 = 1 - b^T G^{-1} b  (G=Gram of rho_{1/k}, b_k=<1,rho_{1/k}>, ||1||^2=1).
# Substitute u=1/x in (1,inf): <f,g> = int_1^inf f(1/u)g(1/u)/u^2 du,  rho_{1/k}(1/u)={u/k}-{u}/k.
import numpy as np
U=4000.0; M=4_000_000            # integrate u in (1,U), fine grid (1/u^2 decay => tail ~1/U)
u=np.linspace(1.0,U,M); du=u[1]-u[0]; w=1.0/u**2
def rho(k):                      # rho_{1/k}(1/u) = {u/k} - {u}/k
    return (u/k - np.floor(u/k)) - (u-np.floor(u))/k
Nmax=30
R=[rho(k) for k in range(1,Nmax+1)]
one_ip=[np.sum(R[k]*w)*du for k in range(Nmax)]      # b_k = <1, rho_{1/k}>
# Gram
G=np.zeros((Nmax,Nmax))
for i in range(Nmax):
    for j in range(i,Nmax):
        G[i,j]=G[j,i]=np.sum(R[i]*R[j]*w)*du
print("N    d_N^2       d_N       (RH: d_N -> 0)   BDLS-predicted ~ (2+gamma-log4pi)/log N ... check trend")
for N in [3,5,8,12,16,20,25,30]:
    GN=G[:N,:N]; bN=np.array(one_ip[:N])
    c=np.linalg.solve(GN+1e-14*np.eye(N), bN)
    d2=1.0 - bN@c
    print(f"{N:3d}  {d2:.6f}   {np.sqrt(max(d2,0)):.6f}")
# optimal coefficients for N=20: do they track Mobius mu(k)/k structure?
N=20; c=np.linalg.solve(G[:N,:N]+1e-14*np.eye(N), np.array(one_ip[:N]))
from sympy import mobius
print("\noptimal c_k (N=20) vs mu(k)/k-ish (Mobius sign pattern):")
print("  k :  c_k      (sign)   mu(k)")
for k in range(1,11):
    print(f"  {k:2d}: {c[k-1]:+.4f}            {int(mobius(k)):+d}")

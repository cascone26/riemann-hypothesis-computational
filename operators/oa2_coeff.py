import numpy as np
from scipy.integrate import quad
# leading form Q0 = (1/2)∫∫|u|log|u| w'w' on (-1,1), sine-deriv basis. Full spectrum + eigvecs.
def R(j,l,u):
    pl=np.pi*l/2; amp=(np.pi*j/2)*pl
    def integ(a,b):
        if abs(a)<1e-30: return (2-u)*np.cos(b)
        return (np.sin(a*2+b)-np.sin(a*u+b))/a
    return amp*0.5*(integ(np.pi*(j-l)/2,pl*u)+integ(np.pi*(j+l)/2,-pl*u))
def k(u):
    au=abs(u); return au*np.log(au) if au>1e-300 else 0.0
def Qelem(j,l):
    v,_=quad(lambda u:k(u)*(R(j,l,u)+R(l,j,u)),0,2,points=[1e-6],limit=200); return v
def oa2(K):
    Q=np.zeros((K,K))
    for j in range(1,K+1):
        for l in range(j,K+1):
            v=0.5*Qelem(j,l); Q[j-1,l-1]=v; Q[l-1,j-1]=v   # factor 1/2 for (1/2)|t|log|t|
    val,vec=np.linalg.eigh(Q)   # mass=I (||w||^2=sum c^2); eigenvectors orthonormal
    jj=np.arange(1,K+1)
    b=np.where(jj%2==1, 4.0/(jj*np.pi), 0.0)   # ∫w_j = 4/(jπ) odd else 0
    c0=vec[:,0]; lam0=val[0]
    cb0=c0@b
    # O(a) coeff check = (7/4)*(∫w0)^2  (||w0||=1)
    oa1=(7.0/4.0)*cb0**2
    # O(a^2) = (7/4)^2 * (c0·b)^2 * sum_{n>=1} (c_n·b)^2/(lam0-lam_n)
    S=0.0
    for n in range(1,K):
        cbn=vec[:,n]@b
        S += cbn**2/(lam0-val[n])
    oa2v=(7.0/4.0)**2 * cb0**2 * S
    return oa1,oa2v,lam0
print("O(a) and O(a^2) coefficients of Suzuki λ_a small-a expansion (2nd-order perturbation):",flush=True)
for K in [40,64,96,128]:
    o1,o2,l0=oa2(K)
    print(f"  K={K:4d}:  O(a)={o1:.5f} (expect ~3.379)   O(a^2)={o2:.5f}",flush=True)

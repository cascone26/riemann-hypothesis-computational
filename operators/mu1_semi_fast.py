import numpy as np
from scipy.integrate import quad
def R(j,l,u):
    pj=np.pi*j/2; pl=np.pi*l/2; amp=pj*pl
    def integ(a,b):
        if abs(a)<1e-30: return (2-u)*np.cos(b)
        return (np.sin(a*2+b)-np.sin(a*u+b))/a
    return amp*0.5*(integ(np.pi*(j-l)/2,pl*u)+integ(np.pi*(j+l)/2,-pl*u))
def k(u):
    au=abs(u); return au*np.log(au) if au>1e-300 else 0.0
def Qelem(j,l):
    f=lambda u: k(u)*(R(j,l,u)+R(l,j,u))
    val,_=quad(f,0,2,points=[1e-6],limit=200); return val
def lam(K):
    Q=np.zeros((K,K))
    for j in range(1,K+1):
        for l in range(j,K+1):
            v=Qelem(j,l); Q[j-1,l-1]=v; Q[l-1,j-1]=v
    return np.linalg.eigvalsh(Q)[0]
print("semi-analytic (M-free, machine-prec elements, scipy quad) K-convergence:",flush=True)
data=[]
for K in [16,24,32,48,64,96,128]:
    lm=lam(K); data.append((K,lm)); print(f"  K={K:4d}  lambda={lm:.9f}  mu1={1+0.5*lm:.9f}",flush=True)
(K2,L2),(K3,L3)=data[-2],data[-1]
for pw in [1.0,2.0]:
    c=(L3-L2)/(K3**-pw - K2**-pw); linf=L3-c*K3**-pw
    print(f"  extrap 1/K^{pw} (last 2): lam_inf={linf:.7f}  mu1={1+0.5*linf:.7f}",flush=True)

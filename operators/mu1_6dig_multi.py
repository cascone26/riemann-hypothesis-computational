import numpy as np, mpmath as mp
from scipy.integrate import quad
from scipy.optimize import curve_fit
mp.mp.dps=30
def R(j,l,u):
    pl=np.pi*l/2; amp=(np.pi*j/2)*pl
    def integ(a,b):
        if abs(a)<1e-30: return (2-u)*np.cos(b)
        return (np.sin(a*2+b)-np.sin(a*u+b))/a
    return amp*0.5*(integ(np.pi*(j-l)/2,pl*u)+integ(np.pi*(j+l)/2,-pl*u))
def k(u):
    au=abs(u); return au*np.log(au) if au>1e-300 else 0.0
def Qel(j,l):
    v,_=quad(lambda u:k(u)*(R(j,l,u)+R(l,j,u)),0,2,points=[1e-6],limit=200); return v
def lam(K):
    Q=np.zeros((K,K))
    for j in range(1,K+1):
        for l in range(j,K+1):
            v=Qel(j,l); Q[j-1,l-1]=v; Q[l-1,j-1]=v
    return np.linalg.eigvalsh(Q)[0]
old={16:-1.492180761,24:-1.505418582,32:-1.512197338,48:-1.519107743,64:-1.522623031,96:-1.526188585,128:-1.527994546,160:-1.529087218,224:-1.530346021,320:-1.531298680}
for K in [200,280,400]:
    old[K]=lam(K); print(f"  K={K} lambda={old[K]:.9f}",flush=True)
K=np.array(sorted(old),float); L=np.array([old[k] for k in sorted(old)])
# multi-term fit L = Linf + c1*K^-p + c2*K^-2p
f=lambda k,Linf,c1,c2,p: Linf + c1*k**(-p) + c2*k**(-2*p)
popt,_=curve_fit(f,K,L,p0=[-1.534,0.5,0.1,0.95],maxfev=500000)
Linf=popt[0]; mu1=1+0.5*Linf
print(f"\nmulti-term fit: Linf={Linf:.8f} p={popt[3]:.4f}  =>  mu1={mu1:.8f}",flush=True)
print(f"  identify(mu1): {mp.identify(mp.mpf(str(round(mu1,7))),['euler','log(2)','log(pi)','pi','catalan','zeta(3)'])}",flush=True)
print(f"  identify(lam): {mp.identify(mp.mpf(str(round(Linf,7))),['euler','log(2)','log(pi)','pi'])}",flush=True)

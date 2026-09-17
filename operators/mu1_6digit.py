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
def Qelem(j,l):
    v,_=quad(lambda u:k(u)*(R(j,l,u)+R(l,j,u)),0,2,points=[1e-6],limit=200); return v
def lam(K):
    Q=np.zeros((K,K))
    for j in range(1,K+1):
        for l in range(j,K+1):
            v=Qelem(j,l); Q[j-1,l-1]=v; Q[l-1,j-1]=v
    return np.linalg.eigvalsh(Q)[0]
# combine old + new K points
old={16:-1.492180761,24:-1.505418582,32:-1.512197338,48:-1.519107743,64:-1.522623031,96:-1.526188585,128:-1.527994546}
Ks=[160,224,320]
for K in Ks:
    old[K]=lam(K); print(f"  K={K} lambda={old[K]:.9f}",flush=True)
K=np.array(sorted(old)); L=np.array([old[k] for k in K])
f=lambda k,Linf,c,p: Linf+c*k**(-p)
popt,_=curve_fit(f,K.astype(float),L,p0=[-1.534,0.57,0.95],maxfev=200000)
Linf=popt[0]; mu1=1+0.5*Linf
print(f"\nfull-fit Linf={Linf:.8f} p={popt[2]:.4f}  =>  mu1={mu1:.8f}",flush=True)
# richardson per adjacent using fitted p
p=popt[2]
rich=[ (L[i+1]-((L[i+1]-L[i])/(K[i+1]**-p-K[i]**-p))*K[i+1]**-p) for i in range(len(K)-1)]
print(f"richardson tail: {[f'{x:.7f}' for x in rich[-4:]]}",flush=True)
print(f"identify(mu1,6dig): {mp.identify(mp.mpf(str(round(mu1,7))),['euler','log(2)','log(pi)','pi','catalan','zeta(3)','sqrt(2)'])}",flush=True)

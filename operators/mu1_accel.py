import numpy as np, mpmath as mp
mp.mp.dps=30
# existing semi-analytic K-sequence (M-free, machine-prec elements)
K=np.array([16,24,32,48,64,96,128],float)
L=np.array([-1.492180761,-1.505418582,-1.512197338,-1.519107743,-1.522623031,-1.526188585,-1.527994546])
# 1) nonlinear fit L = Linf + c*K^{-p}
from scipy.optimize import curve_fit
f=lambda k,Linf,c,p: Linf + c*k**(-p)
popt,_=curve_fit(f,K,L,p0=[-1.53,1.0,1.0],maxfev=100000)
print(f"nonlinear fit L=Linf+c*K^-p:  Linf={popt[0]:.7f}  c={popt[1]:.4f}  p={popt[2]:.4f}",flush=True)
print(f"  => mu1 = {1+0.5*popt[0]:.7f}",flush=True)
# 2) Richardson assuming p (from fit) on consecutive triples
p=popt[2]
def rich2(i):
    (k1,l1),(k2,l2)=(K[i],L[i]),(K[i+1],L[i+1])
    c=(l2-l1)/(k2**-p-k1**-p); return l2-c*k2**-p
print("  Richardson(p) per adjacent pair:",[f"{rich2(i):.6f}" for i in range(len(K)-1)],flush=True)
# 3) Wynn epsilon (Shanks) on the sequence
def shanks(s):
    return [(s[i+2]*s[i]-s[i+1]**2)/(s[i+2]-2*s[i+1]+s[i]) for i in range(len(s)-2)]
s1=shanks(list(L)); s2=shanks(s1) if len(s1)>2 else []
print(f"  Shanks once: {[f'{x:.6f}' for x in s1]}",flush=True)
if s2: print(f"  Shanks twice: {[f'{x:.6f}' for x in s2]}",flush=True)
Linf=popt[0]
print(f"\n  best Linf~{Linf:.6f}  mu1~{1+0.5*Linf:.6f}",flush=True)
print(f"  identify(mu1): {mp.identify(mp.mpf(str(1+0.5*Linf)),['euler','log(2)','log(pi)','pi','catalan'])}",flush=True)

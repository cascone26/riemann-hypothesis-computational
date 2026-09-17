import mpmath as mp, numpy as np
mp.mp.dps=25
def R(j,l,u):
    pj=mp.pi*j/2; pl=mp.pi*l/2; amp=pj*pl
    def integ(a,b):
        if abs(a)<mp.mpf('1e-30'): return (2-u)*mp.cos(b)
        return (mp.sin(a*2+b)-mp.sin(a*u+b))/a
    return amp*mp.mpf('0.5')*(integ(mp.pi*(j-l)/2,pl*u)+integ(mp.pi*(j+l)/2,-pl*u))
def k(u):
    au=abs(u); return au*mp.log(au) if au>0 else mp.mpf(0)
def Qelem(j,l):
    f=lambda u: k(u)*(R(j,l,u)+R(l,j,u))
    return mp.quad(f,[0,mp.mpf('1e-4'),2])
def lam(K):
    Q=mp.matrix(K,K)
    for j in range(1,K+1):
        for l in range(j,K+1):
            v=Qelem(j,l); Q[j-1,l-1]=v; Q[l-1,j-1]=v
    return float(min(mp.eigsy(Q,eigvals_only=True)))
print("semi-analytic (M-free, machine-prec elements) K-convergence:",flush=True)
data=[]
for K in [16,24,32,48,64]:
    lm=lam(K); data.append((K,lm)); print(f"  K={K:3d}  lambda={lm:.8f}  mu1={1+0.5*lm:.8f}",flush=True)
# extrapolate K->inf : lam(K) = lam_inf + c*K^{-p}
K_=np.array([d[0] for d in data],float); L_=np.array([d[1] for d in data])
# fit last 3 to lam_inf + c/K^p via log of successive diffs
d1=np.diff(L_); 
# 1/K extrapolation (Richardson on last 3)
def rich(K1,L1,K2,L2,K3,L3):
    # assume L = Linf + c/K ; solve
    c=(L2-L1)/(1/K2-1/K1); return L1 - c/K1
print(f"\n  1/K Richardson (last 3): lam_inf~{rich(*[v for pair in data[-3:] for v in pair]):.6f}  mu1~{1+0.5*rich(*[v for pair in data[-3:] for v in pair]):.6f}",flush=True)

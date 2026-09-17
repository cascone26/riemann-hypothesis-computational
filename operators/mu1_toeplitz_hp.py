import numpy as np, mpmath as mp
from scipy.linalg import matmul_toeplitz
mp.mp.dps=30
def kernel(u):
    au=np.abs(u); out=np.zeros_like(au); nz=au>1e-300; out[nz]=au[nz]*np.log(au[nz]); return out
def lam_min(M,K):
    x=np.linspace(-1,1,M+1); h=x[1]-x[0]
    d=np.arange(-M,M+1)*h; gv=kernel(d)          # even
    col=gv[M:2*M+1]                                # G[:,0]=f(i)=gv[i+M], i=0..M
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+1),kk)/2.0
    Vp=(kk*np.pi/2.0)*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]                                # (M+1,K)
    GW=matmul_toeplitz((col,col),Wp)               # symmetric Toeplitz @ Wp, FFT O(MK logM)
    Q=Wp.T@GW; Q=(Q+Q.T)/2
    return np.linalg.eigvalsh(Q)[0]
print("FFT-Toeplitz high-precision lambda_min^(1) (|t|log|t| deriv form on (-1,1)):",flush=True)
res=[]
for M,K in [(8000,200),(16000,300),(30000,400),(50000,500)]:
    lm=lam_min(M,K); res.append((M,K,lm))
    print(f"  M={M:6d} K={K:4d}  lambda_min^(1)={lm:.9f}  mu1={1+0.5*lm:.9f}",flush=True)
# Richardson in M (last two, fixed-ish K trend)
lm_best=res[-1][2]; mu1=1+0.5*lm_best
print(f"\nbest lambda_min^(1)={lm_best:.9f}  mu1={mu1:.9f}",flush=True)
# PSLQ inverse-symbolic on mu1 and lambda
print("mpmath.identify(mu1):",mp.identify(mp.mpf(str(mu1)),['euler','log(2)','log(pi)','pi','catalan']),flush=True)
print("mpmath.identify(lambda):",mp.identify(mp.mpf(str(lm_best)),['euler','log(2)','log(pi)','pi']),flush=True)

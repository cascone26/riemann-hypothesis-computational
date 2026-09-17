import numpy as np, mpmath as mp
mp.mp.dps=30
def kernel(u):
    au=np.abs(u); out=np.zeros_like(au); nz=au>1e-300; out[nz]=au[nz]*np.log(au[nz]); return out
def toep_matvec(col, X):   # symmetric Toeplitz (first col `col`, len M1) times X (M1,K), FFT
    M1=len(col); L=1
    while L < 2*M1: L*=2
    c=np.zeros(L); c[:M1]=col; c[L-M1+1:]=col[1:][::-1]
    Cf=np.fft.rfft(c)
    Xp=np.zeros((L,X.shape[1])); Xp[:M1]=X
    Y=np.fft.irfft(np.fft.rfft(Xp,axis=0)*Cf[:,None],n=L,axis=0)
    return Y[:M1]
def lam_min(M,K):
    x=np.linspace(-1,1,M+1); h=x[1]-x[0]
    col=kernel(np.arange(0,M+1)*h)            # G[:,0]=f(i)=kernel(i*h)
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+1),kk)/2.0
    Vp=(kk*np.pi/2.0)*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]
    GW=toep_matvec(col,Wp)
    Q=Wp.T@GW; Q=(Q+Q.T)/2
    return np.linalg.eigvalsh(Q)[0]
print("FFT (hand-rolled) high-precision lambda_min^(1):",flush=True)
res=[]
for M,K in [(8000,250),(20000,350),(50000,450),(100000,550)]:
    lm=lam_min(M,K); res.append((M,K,lm))
    print(f"  M={M:7d} K={K:4d}  lam={lm:.10f}  mu1={1+0.5*lm:.10f}",flush=True)
lm=res[-1][2]; mu1=1+0.5*lm
print(f"\nbest lam={lm:.10f}  mu1={mu1:.10f}",flush=True)
print("identify(mu1):",mp.identify(mp.mpf(str(mu1)),['euler','log(2)','log(pi)']),flush=True)
print("identify(lam):",mp.identify(mp.mpf(str(lm)),['euler','log(2)','log(pi)']),flush=True)

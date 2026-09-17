import numpy as np, mpmath as mp
mp.mp.dps=20
def kernel(u):
    au=np.abs(u); out=np.zeros_like(au); nz=au>1e-300; out[nz]=au[nz]*np.log(au[nz]); return out
def toep_matvec(col,X):
    M1=len(col); L=1
    while L<2*M1: L*=2
    c=np.zeros(L); c[:M1]=col; c[L-M1+1:]=col[1:][::-1]
    Cf=np.fft.rfft(c); Xp=np.zeros((L,X.shape[1])); Xp[:M1]=X
    return np.fft.irfft(np.fft.rfft(Xp,axis=0)*Cf[:,None],n=L,axis=0)[:M1]
def eigfun(M,K):
    x=np.linspace(-1,1,M+1); h=x[1]-x[0]
    col=kernel(np.arange(0,M+1)*h)
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+1),kk)/2.0
    Vp=(kk*np.pi/2.0)*np.cos(arg); w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]; Q=Wp.T@toep_matvec(col,Wp); Q=(Q+Q.T)/2
    val,vec=np.linalg.eigh(Q)
    return val[0],vec[:,0],kk
print("VERIFY O(a) coeff = (7/4)*(∫w0)^2/||w0||^2  [w0 = mu1-eigenfunction]",flush=True)
print("empirical O(a) coeff (from lambda_a small-a fit) ≈ 3.28–3.42",flush=True)
for M,K in [(8000,200),(20000,300),(40000,400)]:
    lam,c,kk=eigfun(M,K)
    # ∫_{-1}^1 w0 = sum_{k odd} c_k * 4/(k pi);  ||w0||^2 = sum c_k^2
    integ=np.sum(c[kk%2==1]*4.0/(kk[kk%2==1]*np.pi))
    norm2=np.sum(c**2)
    ratio=integ**2/norm2
    oa=(7.0/4.0)*ratio
    print(f"  M={M:6d} K={K:4d}: (∫w0)^2/||w0||^2={ratio:.5f}  O(a) coeff=(7/4)*ratio={oa:.4f}",flush=True)

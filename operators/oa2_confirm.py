import numpy as np
from suzuki_screw_lambda_a import screw_g_no_primes, vonmangoldt_upto
def toep_matvec(col,X):
    M1=len(col); L=1
    while L<2*M1: L*=2
    c=np.zeros(L); c[:M1]=col; c[L-M1+1:]=col[1:][::-1]
    Cf=np.fft.rfft(c); Xp=np.zeros((L,X.shape[1])); Xp[:M1]=X
    return np.fft.irfft(np.fft.rfft(Xp,axis=0)*Cf[:,None],n=L,axis=0)[:M1]
def lam_arch(a,M,K):
    Lam=vonmangoldt_upto(4)  # dummy, primes empty for small a
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    col=np.array([screw_g_no_primes(i*h,Lam,4) for i in range(M+1)])
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg); w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]; Q=Wp.T@toep_matvec(col,Wp); Q=(Q+Q.T)/2
    return np.linalg.eigvalsh(Q/a)[0]
A=0.5*(np.log(2*np.pi)+np.euler_gamma-1)
aa=[0.02,0.03,0.04,0.05,0.06,0.08,0.10]
lam=[]
print("high-res arch λ_a (FFT-Toeplitz, M=20000) at small a:",flush=True)
for a in aa:
    l=lam_arch(a,20000,90); lam.append(l); print(f"  a={a:.3f} λ={l:.7f}",flush=True)
a=np.array(aa); lam=np.array(lam); res=lam-np.log(1/a)
for deg in [2,3]:
    c=np.polyfit(a,res,deg)
    print(f"  {deg+1}-term fit: O(a^2)={c[-3]:.4f}  O(a)={c[-2]:.4f}  mu1={c[-1]+2*A+1:.5f}",flush=True)
print("  predicted: O(a^2)=-0.196, O(a)=3.379, mu1=0.2331",flush=True)

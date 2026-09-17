import numpy as np
A=0.5*(np.log(2*np.pi)+np.euler_gamma-1)
def k_lead(t):
    at=abs(t)
    return (0.5*at*np.log(at) + A*at) if at>1e-300 else 0.0
def toep_matvec(col,X):
    M1=len(col); L=1
    while L<2*M1: L*=2
    c=np.zeros(L); c[:M1]=col; c[L-M1+1:]=col[1:][::-1]
    Cf=np.fft.rfft(c); Xp=np.zeros((L,X.shape[1])); Xp[:M1]=X
    return np.fft.irfft(np.fft.rfft(Xp,axis=0)*Cf[:,None],n=L,axis=0)[:M1]
def lam_lead(a,M,K):
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    col=np.array([k_lead(i*h) for i in range(M+1)])
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg); w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]; Q=Wp.T@toep_matvec(col,Wp); Q=(Q+Q.T)/2
    return np.linalg.eigvalsh(Q/a)[0]
aa=[0.02,0.03,0.04,0.05,0.06,0.08,0.10]
lam=np.array([lam_lead(a,20000,90) for a in aa]); a=np.array(aa)
res=lam-np.log(1/a)
c=np.polyfit(a,res,3)
print(f"LEADING-only kernel ½|t|log|t|+A|t|:  fit const C={c[-1]:.5f}  O(a)={c[-2]:.4f}",flush=True)
print(f"  => mu1_lead = C+2A+1 = {c[-1]+2*A+1:.5f}",flush=True)
print(f"  compare: full-arch fit μ₁=0.237 ; unit-interval semi/grid μ₁=0.233",flush=True)
# also: directly, the leading form should give lambda_a = log(1/a) + (1/2)lam_min1 - 2A exactly.
# so lam_min1_implied = 2*(C+2A) ... res->C = (1/2)lam_min1 - 2A  => lam_min1 = 2(C+2A)
print(f"  => implied lambda_min^(1) = 2(C+2A) = {2*(c[-1]+2*A):.5f}  (unit-interval semi gave -1.5337)",flush=True)

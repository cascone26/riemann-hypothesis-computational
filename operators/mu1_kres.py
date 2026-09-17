import numpy as np
A=0.5*(np.log(2*np.pi)+np.euler_gamma-1)
def k_lead(t):
    at=abs(t); return (0.5*at*np.log(at)+A*at) if at>1e-300 else 0.0
def toep(col,X):
    M1=len(col); L=1
    while L<2*M1: L*=2
    c=np.zeros(L); c[:M1]=col; c[L-M1+1:]=col[1:][::-1]
    Cf=np.fft.rfft(c); Xp=np.zeros((L,X.shape[1])); Xp[:M1]=X
    return np.fft.irfft(np.fft.rfft(Xp,axis=0)*Cf[:,None],n=L,axis=0)[:M1]
def lam(a,M,K):
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    col=np.array([k_lead(i*h) for i in range(M+1)])
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg); w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]; Q=Wp.T@toep(col,Wp); Q=(Q+Q.T)/2
    return np.linalg.eigvalsh(Q/a)[0]
# leading kernel: lambda_a = log(1/a)+const EXACTLY, so const from ONE a is enough. Use a=0.05.
print("leading-kernel μ₁ vs basis K (at a=0.05, M=20000) — does it converge to 0.233?",flush=True)
for K in [90,140,200,280]:
    l=lam(0.05,20000,K); const=l-np.log(1/0.05); mu1=const+2*A+1
    print(f"  K={K:4d}: λ_a(0.05)={l:.6f}  μ₁={mu1:.5f}  λ_min⁽¹⁾={2*(const+2*A):.5f}",flush=True)
print("  unit-interval semi K→∞: μ₁=0.2331 (λ_min⁽¹⁾=-1.5337)",flush=True)

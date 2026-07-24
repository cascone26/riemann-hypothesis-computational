# Sign-definite modular parent = Polya Phi(u): Xi(z)=int Phi(u)e^{izu}du, Phi>0 & even (theta FE).
import numpy as np
def Phi(u, nmax=40):
    u=np.atleast_1d(np.asarray(u,dtype=float)); n=np.arange(1,nmax+1)[:,None]
    e4=np.exp(4*u)[None,:]
    T=(2*np.pi**2*n**4*np.exp(9*u)[None,:]-3*np.pi*n**2*np.exp(5*u)[None,:])*np.exp(-np.pi*n**2*e4)
    return T.sum(0)
print("(1) Phi positive & even (Phi>0 is the theta-FE 'miracle'):")
for u in [-1.0,-0.5,-0.2,0.0,0.2,0.5,1.0]:
    print(f"   Phi({u:+.1f})={Phi(u)[0]:+.6e}   Phi(-u)-Phi(u)={ (Phi(-u)-Phi(u))[0]:+.2e}")
# cosine transform F(z)=2 int_0^U e^{lam u^2}Phi(u)cos(zu)du (Phi even => real Xi)
U=3.0; M=6000; u=np.linspace(0,U,M+1); du=u[1]-u[0]; w=np.full(M+1,du); w[0]=w[-1]=du/2
phi=Phi(u)
def Xitr(z,lam=0.0):
    return 2*np.sum(w*np.exp(lam*u*u)*phi*np.cos(z*u))
F0=Xitr(0.0)
print(f"\n(2) cosine transform of Phi = Xi ? vanish at zeros gamma_n:")
for z,lab in [(0.0,"0"),(17.5,"17.5 off"),(14.134725,"g1"),(21.02204,"g2"),(25.010858,"g3")]:
    print(f"   F({lab:9s})={Xitr(z):+.6e}   |F|/|F0|={abs(Xitr(z)/F0):.2e}")
print("\n(3) de Bruijn-Newman heat flow: real zero near g1 vs lambda (RH<=>Lambda<=0; Rodgers-Tao: Lambda>=0):")
def find_zero(lam, lo=13.9, hi=14.4):
    flo=Xitr(lo,lam)
    for _ in range(60):
        mid=(lo+hi)/2; fm=Xitr(mid,lam)
        if (fm>0)==(flo>0): lo,flo=mid,fm
        else: hi=mid
    return (lo+hi)/2
for lam in [-0.05,0.0,0.05,0.1]:
    print(f"   lambda={lam:+.2f}: zero near g1 at z={find_zero(lam):.6f}")
print("Lambda=0 critical; zeros stay real for lambda>=Lambda. RH<=>Lambda<=0, Lambda>=0 PROVEN (Rodgers-Tao 2018).")

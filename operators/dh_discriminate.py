import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
import mpmath as mp
mp.mp.dps=18
import numpy as np
from flint import arb, ctx
from _dh_negative_control import lambda_D

def psiR(x,L,nu,cond):
    x=mp.mpf(x); nu=mp.mpf(nu); L=mp.mpf(L); logc=mp.log(mp.pi/cond)
    hp=lambda r: mp.re(mp.digamma(nu+1j*r/2))-logc
    def Sr(r):
        B=2*mp.pi*x/L; den=B*B-r*r
        if abs(den)<mp.mpf(10)**-18: r=r+mp.mpf(10)**-12; den=B*B-r*r
        return (mp.cos(r*L)-mp.cos(2*mp.pi*x))*B/den
    f=lambda r: hp(r)*Sr(r); B=abs(2*mp.pi*x/L) if x!=0 else mp.mpf(1)
    per=mp.pi/L; R=mp.mpf(40); pts=[mp.mpf(0)]; k=2
    while k*per<R: pts.append(k*per); k+=2
    pts.append(R); allp=sorted(set([float(p) for p in pts]+[float(B)]))
    val=2*mp.quad(f,allp)+2*B*mp.cos(2*mp.pi*x)*mp.quad(lambda r: hp(r)/r**2,[R,4*R,mp.inf])
    return -val/(2*mp.pi**2)

def arch(N,L,nu,cond):
    v={k:psiR(k,L,nu,cond) for k in range(N+1)}; h=mp.mpf(10)**-5
    dg={n:(psiR(n+h,L,nu,cond)-psiR(n-h,L,nu,cond))/(2*h) for n in range(N+1)}
    A=np.zeros((2*N+1,2*N+1))
    for i in range(2*N+1):
        n=i-N
        for j in range(2*N+1):
            m=j-N
            A[i,j]=float(dg[abs(n)]) if n==m else float((v[abs(m)]-v[abs(n)])/(m-n))
    return A

def prime(N,L,w,pos):
    A=np.zeros((2*N+1,2*N+1))
    for i in range(2*N+1):
        n=i-N
        for j in range(2*N+1):
            m=j-N; s=0.0
            for wi,y in zip(w,pos):
                q=2*(1-y/L)*math.cos(2*math.pi*n*y/L) if n==m else (math.sin(2*math.pi*m*y/L)-math.sin(2*math.pi*n*y/L))/(math.pi*(n-m))
                s+=wi*q
            A[i,j]=s
    return A

c=14;N=8;L=float(mp.log(c))
def vm(N):
    lam=[0.0]*(N+1)
    for p in range(2,N+1):
        if lam[p]==0:
            lp=math.log(p);pk=p
            while pk<=N: lam[pk]=lp;pk*=p
    return lam
VM=vm(c); zpos=[math.log(n) for n in range(2,c+1) if VM[n]>0]; zw=[VM[n]/math.sqrt(n) for n in range(2,c+1) if VM[n]>0]
LD=lambda_D(c); dpos=[math.log(n) for n in range(2,c+1)]; dw=[LD[n]/math.sqrt(n) for n in range(2,c+1)]
ctx.prec=200; La=arb(c).log(); PI=arb.pi(); sp2=16*PI*PI; l2=La*La; pref=32*La*(La/4).sinh()**2
W02=np.zeros((2*N+1,2*N+1))
for i in range(2*N+1):
    n=i-N
    for j in range(2*N+1):
        m=j-N; W02[i,j]=float((pref*(l2-sp2*m*n)/((l2+sp2*m*m)*(l2+sp2*n*n))).mid())
print(f"c={c} N={N}: zeta vs D-H Weil min-eig (validated-numeric arch)",flush=True)
tauZ=W02+arch(N,L,mp.mpf(1)/4,1)-prime(N,L,zw,zpos)
tauD=arch(N,L,mp.mpf(3)/4,5)-prime(N,L,dw,dpos)   # NO pole (entire)
eZ=np.linalg.eigvalsh(0.5*(tauZ+tauZ.T)); eD=np.linalg.eigvalsh(0.5*(tauD+tauD.T))
print(f"  zeta: min-eig={eZ[0]:+.5e} n_neg={int((eZ<-1e-9).sum())}",flush=True)
print(f"  D-H : min-eig={eD[0]:+.5e} n_neg={int((eD<-1e-9).sum())}",flush=True)
print(f"  => {'DISCRIMINATED (D-H<zeta)' if eD[0]<eZ[0]-1e-6 else 'inconclusive'}",flush=True)

import numpy as np, mpmath as mp, math
mp.mp.dps=25; PI=np.pi
atil=np.load('E11_an.npy'); NMc=len(atil)-1
Zc=np.load('/tmp/E11_zeros.npy')
# Lambda_E(n) via log-derivative: atil_n log n = sum_{d|n} Lambda_E(d) atil_{n/d}, atil_1=1
NM=2000
at=[0.0]+[float(atil[n]) for n in range(1,NM+1)]
Lam=[0.0]*(NM+1)
for n in range(2,NM+1):
    tot=at[n]*math.log(n); d=1
    while d*d<=n:
        if n%d==0:
            if 1<d<n: tot-=Lam[d]*at[n//d]
            e=n//d
            if 1<e<n and e!=d: tot-=Lam[e]*at[n//e]
        d+=1
    Lam[n]=tot
# Euler-product check: Lambda_E ~0 on composites?
import sympy
comps=[round(Lam[n],4) for n in range(2,40) if len(sympy.factorint(n))>1]
print(f"Lambda_E on small composites (Euler product => ~0): max|.|={max(abs(x) for x in comps):.2e}",flush=True)
# arch: Omega_E11(r)=Re psi(1+ir)-log(2pi)+0.5 log(11)
R=200.0;NG=40000;rg=np.linspace(-R,R,NG)
OM=np.array([float(mp.re(mp.digamma(mp.mpc(1,r))))-np.log(2*PI)+0.5*np.log(11) for r in rg])
def ph(a,r): s=np.sqrt(2*PI); return s/2*(np.exp(-(r-a)**2/2)+np.exp(-(r+a)**2/2))
def arch(i,j): return np.trapezoid(ph(i,rg)*ph(j,rg)*OM,rg)/(2*PI)
def g(i,j,y):
    s=i+j;d=i-j
    return np.exp(-y**2/4)*(np.sqrt(PI)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
NS=np.arange(2,NM+1);LN=np.log(NS);W=np.array(Lam[2:])/np.sqrt(NS)
def prime(i,j): return np.sum(W*g(i,j,LN))
def Qa(i,j): return 2*(arch(i,j)-prime(i,j))
def Qg(i,j): return 2*np.sum(ph(i,Zc)*ph(j,Zc))
A=[10,14,18,22,26]
Qam=np.array([[Qa(A[i],A[j]) for j in range(len(A))] for i in range(len(A))])
Qgm=np.array([[Qg(A[i],A[j]) for j in range(len(A))] for i in range(len(A))])
print(f"E11a genuine degree-2 Weil form (5x5, a={A}), {len(Zc)} zeros:")
print(f"  diag Q_arith = {np.round(np.diag(Qam),4).tolist()}")
print(f"  diag Q_gram  = {np.round(np.diag(Qgm),4).tolist()}")
print(f"  max|Q_arith-Q_gram| = {np.max(np.abs(Qam-Qgm)):.3e}")
print(f"  min-eig(Q_arith)={np.linalg.eigvalsh(Qam).min():+.4f}  min-eig(Q_gram)={np.linalg.eigvalsh(Qgm).min():+.4f}  (GRH=>PD, n_neg=0)")

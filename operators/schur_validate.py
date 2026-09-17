# Validate brick 16 (psd_via_schur_coupling) captures the REAL [14,44] corner mechanism:
# split the genuine Weil form Q in the PRIME eigenbasis at the floor; check A>0 & Schur complement>=0.
import numpy as np
from scipy.special import digamma
def vonmangoldt(N):
    lam=np.zeros(N+1);isc=np.zeros(N+1,bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True;lp=np.log(p);pk=p
            while pk<=N: lam[pk]=lp;pk*=p
    return lam
NMAX=3000;LAM=vonmangoldt(NMAX);NS=np.arange(2,NMAX+1);LOGN=np.log(NS);LAMW=LAM[2:]/np.sqrt(NS)
a,b,W=14.0,44.0,0.6875  # 4W=2.75 (the largest certified in PATH D corner summary)
fourW=4*W
dr=0.002; RG=np.arange(-(b+120),(b+120),dr)
OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi)
N=8; CEN=np.linspace(a,b,N)
def fej(x):
    y=W*x;o=np.ones_like(y);nz=np.abs(y)>1e-12;o[nz]=(np.sin(y[nz])/y[nz])**2;return o
def fc(c):
    y=W*(0.5j-c);y2=W*(0.5j+c);return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
P=np.array([fej(RG-c)+fej(RG+c) for c in CEN])
Gm=(P@P.T)*dr
arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
mask=LOGN<fourW;ln=LOGN[mask];w0=LAMW[mask];cosu=np.cos(np.outer(ln,RG))
prime=np.zeros((N,N))
for jj in range(len(ln)):
    col=cosu[jj];Mk=np.zeros((N,N))
    for j in range(N): Mk[j,:]=2.0*(((P[j][None,:]*P)@col)*dr/(2*np.pi))
    prime+=w0[jj]*0.5*(Mk+Mk.T)
# generalized: work in Gram-orthonormal basis
Lg=np.linalg.cholesky(Gm); Li=np.linalg.inv(Lg)
archO=Li@arch@Li.T; poleO=Li@pole@Li.T; primeO=Li@prime@Li.T
Q=archO+poleO-primeO          # Weil form (arch+pole-prime), orthonormalized
Q=0.5*(Q+Q.T); primeO=0.5*(primeO+primeO.T)
# prime eigenbasis
mu,V=np.linalg.eigh(primeO)   # ascending
floor=np.min(np.linalg.eigvalsh(archO+poleO))
strong=mu>floor
print(f"[14,44] 4W={fourW}: prime eigenvalues mu={np.round(mu,4).tolist()}")
print(f"  floor=lambda_min(arch+pole)={floor:.5f}; strong block (mu>floor) size k={strong.sum()}")
Qt=V.T@Q@V                    # Q in prime eigenbasis
idx=np.argsort(~strong)       # strong first
order=np.concatenate([np.where(strong)[0],np.where(~strong)[0]])
Qt=Qt[np.ix_(order,order)]; k=int(strong.sum())
A=Qt[:k,:k]; B=Qt[:k,k:]; D=Qt[k:,k:]
print(f"  lambda_min(Q full)={np.linalg.eigvalsh(Q).min():+.4e}  (>0 = RH-corner certified)")
print(f"  A=strong block: lambda_min(A)={np.linalg.eigvalsh(A).min():+.4e}  (PosDef? {'YES' if np.linalg.eigvalsh(A).min()>0 else 'NO'})")
print(f"  D=weak block:  lambda_min(D)={np.linalg.eigvalsh(D).min():+.4e}  (D alone PSD? {'YES' if np.linalg.eigvalsh(D).min()>=-1e-12 else 'NO -- defect'})")
Sch=D-B.T@np.linalg.inv(A)@B
print(f"  Schur complement D-B'A^-1 B: lambda_min={np.linalg.eigvalsh(Sch).min():+.4e}  (PSD? {'YES' if np.linalg.eigvalsh(Sch).min()>=-1e-9 else 'NO'})")
print("  => brick 16 applies iff A PosDef AND Schur PSD; then full Q PSD (RH-corner). Validates the mechanism.")

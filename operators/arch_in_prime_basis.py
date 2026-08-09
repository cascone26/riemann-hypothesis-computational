#!/usr/bin/env python3
"""§879: the joint structure — arch in the prime eigenbasis at the corner. Q = arch - prime. In the prime
eigenbasis {w_j} (prime diagonal = mu_j), Q_ij = <w_i,arch w_j> - mu_j delta_ij. Check: (1) does arch(w_j)>=mu_j
per-diagonal for ALL j (not just top)? (2) are off-diagonal arch cross-terms small (=> near-diagonal, bound easy)
or load-bearing? (3) what makes lmin(Q)>0."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def build(a,b,W,N,dr=0.004):
    fourW=4*W; RGRID=np.arange(-(b+120),(b+120),dr); OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
    CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr; arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp)); AP=pole+arch
    gw,gv=eigh(G); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    return B@AP@B, B@prime@B
a,b,W,N=14,44,1.0,13   # corner 4W=4
APw,primew=build(a,b,W,N)
mu,Wv=eigh(primew)  # ascending mu, Wv columns = prime eigenvectors (orthonormal)
# arch in prime eigenbasis
Ap=Wv.T@APw@Wv
Q=APw-primew; margin=eigh(Q,eigvals_only=True)[0]
print(f"§879 [14,44] 4W=4 N={N}: margin={margin:+.3e}")
print("prime eigenvalues mu_j (desc):", np.round(mu[::-1],4))
print("arch(w_j)=Ap_jj (same order desc):", np.round(np.diag(Ap)[::-1],4))
diff=np.diag(Ap)[::-1]-mu[::-1]
print("arch(w_j)-mu_j (desc):", np.round(diff,4), " all>=0?", bool((diff>=-1e-9).all()))
# off-diagonal magnitude vs diagonal
offmax=np.max(np.abs(Ap-np.diag(np.diag(Ap))))
print(f"max |off-diag arch in prime basis| = {offmax:.4f}  vs min diag Q = {np.min(np.diag(Ap)-mu):.4f}")
# Q in prime basis: diag = arch(w_j)-mu_j, offdiag = arch cross terms
Qp=Ap-np.diag(mu)
print("Q diag (arch(w_j)-mu_j) desc:", np.round(np.diag(Qp)[::-1],4))
print(f"Are off-diagonals load-bearing? If |offdiag| ~ diag, yes. offmax={offmax:.4f}, min|diag|={np.min(np.abs(np.diag(Qp))):.4f}")

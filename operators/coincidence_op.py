#!/usr/bin/env python3
"""§882: characterize the coincidence operator M = arch|_S ~= prime|_S on the critical subspace, and confirm the
precise RH-corner condition is Schur(Q_SS) >= 0 (Q_SS itself is only necessary). Is M a recognizable operator
(concentration/prolate)? Compare arch|_S, prime|_S, and the concentration operator A=2 int_window F F restricted to S."""
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
    win=(RGRID>=a)&(RGRID<=b); Aconc=(P[:,win]@P[:,win].T)*dr*2.0  # concentration operator
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp)); AP=pole+arch
    gw,gv=eigh(G); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    return B@AP@B, B@prime@B, B@Aconc@B
for (a,b,W,N,lbl) in [(14,44,1.0,13,'4W=4'),(14,44,0.75,10,'4W=3')]:
    APw,primew,Aconc=build(a,b,W,N); floor=eigh(APw,eigvals_only=True)[0]
    mu,Wv=eigh(primew); idx=np.argsort(mu)[::-1]; mus=mu[idx]; Wv=Wv[:,idx]
    k=int((mus>floor).sum()); Ws=Wv[:,:k]
    archS=Ws.T@APw@Ws; primeS=Ws.T@primew@Ws; concS=Ws.T@Aconc@Ws
    QSS=archS-primeS
    # Schur: full Q, Schur complement of the S block wrt R
    Qw=APw-primew; Qp=Wv.T@Qw@Wv; QSSf=Qp[:k,:k]; QSR=Qp[:k,k:]; QRR=Qp[k:,k:]
    schur=QSSf-QSR@np.linalg.inv(QRR)@QSR.T
    print(f"§882 [{a},{b}] {lbl}: k={k}, floor={floor:.4f}")
    print(f"  arch|_S eigenvalues: {np.round(eigh(archS,eigvals_only=True),4)}")
    print(f"  prime|_S eigenvalues (mu_j): {np.round(np.sort(mus[:k]),4)}")
    print(f"  concentration|_S eigenvalues: {np.round(eigh(concS,eigvals_only=True),4)}")
    print(f"  ||arch|_S - prime|_S||_F = {np.linalg.norm(QSS):.4e}  (coincidence)")
    print(f"  ||arch|_S - concentration|_S||_F = {np.linalg.norm(archS-concS):.4e}  (is M the concentration op?)")
    print(f"  Q_SS min eig = {eigh(QSS,eigvals_only=True)[0]:+.4e};  SCHUR min eig = {eigh(schur,eigvals_only=True)[0]:+.4e} (= margin, the TIGHT condition)")
    print()
print("RH-corner precise: Schur(Q_SS) >= 0 (margin). Q_SS>=0 necessary. Is M=arch|_S=prime|_S the concentration operator?")

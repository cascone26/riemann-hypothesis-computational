#!/usr/bin/env python3
"""§880: the corner Schur complement. Split (prime eigenbasis) into strong block S (prime eig > floor) and robust
R (rest). Q PSD <=> Schur = Q_SS - Q_SR Q_RR^-1 Q_RS >= 0. Compute it, its eigenvalues, and what makes it >=0."""
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
for (a,b,W,N,lbl) in [(14,44,1.0,13,'4W=4'),(14,44,0.75,10,'4W=3')]:
    APw,primew=build(a,b,W,N); mu,Wv=eigh(primew)
    floor=eigh(APw,eigvals_only=True)[0]
    Q=APw-primew; Qp=Wv.T@Q@Wv  # Q in prime eigenbasis
    idx=np.argsort(mu)[::-1]     # desc
    Qp=Qp[np.ix_(idx,idx)]; mus=mu[idx]
    k=int((mus>floor).sum())     # strong block size
    S=slice(0,k); R=slice(k,N)
    QSS=Qp[S,S]; QSR=Qp[S,R]; QRR=Qp[R,R]
    schur=QSS - QSR@np.linalg.inv(QRR)@QSR.T
    print(f"§880 [{a},{b}] {lbl} N={N}: floor={floor:.4f}, strong block k={k} (prime eig>floor)")
    print(f"  margin(Q)={eigh(Q,eigvals_only=True)[0]:+.3e}")
    print(f"  Q_SS eigenvalues: {np.round(eigh(QSS,eigvals_only=True),4)}  (min {eigh(QSS,eigvals_only=True)[0]:+.4f})")
    print(f"  Q_RR min eig: {eigh(QRR,eigvals_only=True)[0]:+.4f} (robust block, should be >0)")
    print(f"  SCHUR eigenvalues: {np.round(eigh(schur,eigvals_only=True),6)}")
    print(f"  Schur min = {eigh(schur,eigvals_only=True)[0]:+.3e}  (= margin? this is the RH-hard quantity)")
    print()
print("Q_SS indefinite (near-null diag + coupling) but Schur>=0 => the coupling to R (Q_SR Q_RR^-1 Q_RS) is what")
print("lifts Q_SS to PSD. The Schur min eig is the concrete RH-marginal quantity.")

#!/usr/bin/env python3
"""§885: CHARACTERIZE THE RAZOR. RH-corner ⟺ Schur(Q_SS) = Q_SS − C ⪰ 0, C = Q_SR Q_RR^{-1} Q_RS (the S–R
coupling correction, PSD). The workflow surfaced C as the razor (eats ~97% of Q_SS). Pin down what C IS:
its rank, spectrum, and how much of it is carried by the rank-2 POLE vs the defect D=arch−prime. If C is
low-rank / pole-dominated → the razor is a concrete finite rank-2 obstruction, a big structural simplification."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
np.set_printoptions(precision=3,suppress=True,linewidth=160)

def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)

def blocks(a,b,fourW,N,dr=0.004):
    W=fourW/4; RGRID=np.arange(-(b+120),(b+120),dr)
    OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw,gv=eigh(G); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    archw=B@arch@B; primew=B@prime@B; polew=B@pole@B
    return archw,primew,polew,B

for (a,b,fourW,N) in [(14,44,4.0,17),(14,44,4.0,13)]:
    archw,primew,polew,B=blocks(a,b,fourW,N)
    Q=archw+polew-primew; D=archw-primew
    floor=eigh(archw+polew,eigvals_only=True)[0]
    mu,Wv=eigh(primew); idx=np.argsort(mu)[::-1]; k=int((mu[idx]>floor).sum())
    Vall=Wv[:,idx]  # prime eigenbasis, ordered strong->weak
    Qb=Vall.T@Q@Vall; Db=Vall.T@D@Vall; poleb=Vall.T@polew@Vall  # Q in prime eigenbasis
    S=slice(0,k); R=slice(k,N)
    QSS=Qb[S,S]; QSR=Qb[S,R]; QRR=Qb[R,R]; QRS=Qb[R,S]
    C=QSR@np.linalg.solve(QRR,QRS)  # correction (PSD)
    schur=QSS-C
    mQSS=eigh(QSS,eigvals_only=True)[0]; mSch=eigh(schur,eigvals_only=True)[0]
    Cev=eigh(C,eigvals_only=True)[::-1]
    # how much of C is pole-carried? split coupling into pole part and D part
    poleSR=poleb[S,R]; DSR=Db[S,R]
    C_pole=poleSR@np.linalg.solve(QRR,poleSR.T); C_D=DSR@np.linalg.solve(QRR,DSR.T)
    print(f"\n[{a},{b}] 4W={fourW} N={N}: k={k} (dim S), dim R={N-k}, floor={floor:.4f}")
    print(f"  min Q_SS   = {mQSS:+.3e}   min Schur = {mSch:+.3e}   (razor = full margin)")
    print(f"  C=Q_SR Q_RR^-1 Q_RS eigenvalues (desc): {Cev[:min(6,N-k)]}")
    print(f"  rank(C) numeric (>1e-9): {(Cev>1e-9).sum()}   trace(C)={C.trace():.3e}   trace(Q_SS)={QSS.trace():.3e}")
    print(f"  ||C_pole||={np.linalg.norm(C_pole):.3e}  ||C_from_D||={np.linalg.norm(C_D):.3e}  ||QSR pole part||/||QSR||={np.linalg.norm(poleSR)/np.linalg.norm(QSR):.3f}")
    print(f"  pole rank={np.sum(eigh(polew,eigvals_only=True)>1e-9)}  ||pole||={np.linalg.norm(polew):.3e}  ||D||={np.linalg.norm(D):.3e}")

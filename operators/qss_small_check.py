#!/usr/bin/env python3
"""§883 pre-check: find a SMALL-N corner where Q_SS min eig > arb radii (~1e-6), feasible for arb cert."""
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
ZT=np.array([14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073])
def probe(a,b,W,N,dr=0.004):
    fourW=4*W; RGRID=np.arange(-(b+100),(b+100),dr); OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
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
    APw=B@AP@B; primew=B@prime@B; floor=eigh(APw,eigvals_only=True)[0]
    mu,Wv=eigh(primew); idx=np.argsort(mu)[::-1]; k=int((mu[idx]>floor).sum()); Ws=Wv[:,idx[:k]]
    QSS=Ws.T@(APw-primew)@Ws
    margin=eigh(APw-primew,eigvals_only=True)[0]
    nz=int(((ZT>=a)&(ZT<=b)).sum())
    return k,nz,(eigh(QSS,eigvals_only=True)[0] if k>0 else float('nan')),margin
print(f"{'window':>10} {'4W':>4} {'N':>3} {'#z':>3} {'k(strong)':>9} {'Q_SS min':>11} {'margin':>11}")
for (a,b,fourW,N) in [(14,30,4.0,6),(14,30,4.0,7),(14,30,3.0,5),(14,34,4.0,7),(14,30,5.0,7),(14,26,4.0,6)]:
    W=fourW/4; k,nz,qm,mg=probe(a,b,W,N)
    print(f"{f'[{a},{b}]':>10} {fourW:>4.1f} {N:>3} {nz:>3} {k:>9} {qm:>+11.3e} {mg:>+11.3e}")
print("want: Q_SS min > 2e-6 (certifiable), small N (<=7), and it's a corner (margin tiny, k>=1)")

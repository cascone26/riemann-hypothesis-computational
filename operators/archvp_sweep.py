#!/usr/bin/env python3
"""§877: broad check — does arch(v_p) >= pnorm hold everywhere (v_p=prime top eigvec)? And is arch(v_p)=pnorm
exactly the RH-marginal (margin->0) locus? Sweep many (window,4W,N)."""
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
def probe(a,b,W,N,dr=0.005):
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
    APw=B@AP@B; primew=B@prime@B
    wp,Vp=eigh(primew); pnorm=wp[-1]; vp=Vp[:,-1]; arch_vp=float(vp@APw@vp)
    margin=float(eigh(APw-primew,eigvals_only=True)[0])
    return arch_vp,pnorm,margin
print("§877 broad check arch(v_p)>=pnorm and its relation to margin->0")
fails=0; n=0; mind=1e9
for a in [14,30,60,120,300]:
    for L in [20,30,40]:
        for fourW in [2.0,2.5,3.0,3.5,4.0]:
            W=fourW/4; b=a+L; Ncap=L*2*W/np.pi
            for frac in [0.5,0.7,0.9]:
                N=max(4,round(frac*Ncap))
                if N>22: continue
                avp,pn,mg=probe(a,b,W,N); d=avp-pn; n+=1; mind=min(mind,d)
                if d<-1e-6: fails+=1; print(f"  FAIL [{a},{b}] 4W={fourW} N={N}: arch(vp)-pnorm={d:.4e}, margin={mg:.2e}")
print(f"checked {n} configs. arch(v_p)>=pnorm failures (< -1e-6): {fails}. min(arch(vp)-pnorm)={mind:.4e}")
print("if 0 failures and min~0 at the razor => arch(v_p)>=pnorm holds universally, =0 exactly at RH-margin.")

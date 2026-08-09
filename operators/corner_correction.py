#!/usr/bin/env python3
"""§876: the corner misalignment correction. margin = floor - pnorm + Delta (Delta>0 = misalignment lift).
At the corner floor<pnorm; margin>0 needs Delta >= pnorm-floor. Decompose: v_p = prime top eigvec, v_a = arch
bottom eigvec. Compute arch(v_p) [arch in prime's top direction] - if arch(v_p) >= pnorm then Q(v_p)>=0. Look at
what carries Delta."""
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
def probe(a,b,W,N,dr=0.004):
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
    # generalized eigvecs in G-metric: whiten
    gw,gv=eigh(G); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    APw=B@AP@B; primew=B@prime@B; Qw=APw-primew
    floor,_=eigh(APw); floor=floor[0]
    wp,Vp=eigh(primew); pnorm=wp[-1]; vp=Vp[:,-1]              # prime top (whitened, unit)
    wa,Va=eigh(APw); vfa=Va[:,0]                                # arch bottom
    margin=float(eigh(Qw,eigvals_only=True)[0])
    arch_vp=float(vp@APw@vp)          # arch+pole in prime's top direction
    prime_va=float(vfa@primew@vfa)    # prime in arch's bottom direction
    Q_vp=float(vp@Qw@vp)              # Q on prime-top vector
    cos=abs(float(vp@vfa))
    return floor,pnorm,margin,arch_vp,prime_va,Q_vp,cos
print("="*104)
print("§876 corner correction: does arch(v_p) [arch in prime-top dir] >= pnorm? (=> Q(v_p)>=0). Delta=margin-(floor-pnorm)")
print("="*104)
print(f"{'window':>10} {'4W':>4} {'floor':>7} {'pnorm':>7} {'margin':>11} {'arch(vp)':>9} {'Q(vp)':>9} {'prime(va)':>9} {'Delta':>8} {'cos':>6}")
for (a,b,fourW,N) in [(14,44,2.0,6),(14,44,3.0,10),(14,44,4.0,13),(50,80,4.0,12),(300,330,5.0,14)]:
    W=fourW/4; fl,pn,mg,avp,pva,qvp,cos=probe(a,b,W,N); delta=mg-(fl-pn)
    print(f"{f'[{a},{b}]':>10} {fourW:>4.1f} {fl:>7.4f} {pn:>7.4f} {mg:>+11.3e} {avp:>9.4f} {qvp:>+9.4f} {pva:>9.4f} {delta:>+8.4f} {cos:>6.3f}")
print("="*104)
print("If arch(v_p) >= pnorm consistently => the prime's strongest direction sits where arch is big enough => Q(v_p)>=0.")
print("Delta = the misalignment lift (= pnorm-floor at the corner where margin~0).")

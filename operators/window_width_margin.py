#!/usr/bin/env python3
"""§859: §858 follow-up. The near-null mode vanishes at in-window zeros, positivity from out-window tails.
Predict: WIDER window (more in-window zeros) -> mode must vanish at more points -> residual from further-out
zeros (smaller tails) -> margin THINS. Fixed base a=14, vary width L; fixed 4W; N=0.75*capacity. Track
margin, #in-window zeros, and the nearest out-of-window zero above the window."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.special import digamma
from scipy.linalg import eigh
import mpmath as mp
def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
mp.mp.dps=15; ZT=[]; k=1
while True:
    g=float(mp.im(mp.zetazero(k))); ZT.append(g)
    if g>200: break
    k+=1
ZT=np.array(ZT)
def margin(a,b,W,N,dr=0.004):
    fourW=4*W; RGRID=np.arange(-(b+120),(b+120),dr)
    OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr
    gw=np.linalg.eigvalsh(G); cond=float(gw[-1]/max(gw[0],1e-300))
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    return float(eigh(pole+arch-prime,G,eigvals_only=True)[0]), cond
print("="*88); print("§859 margin vs window WIDTH (base a=14, 4W=4, N=0.75*cap) — does wider window thin the margin?"); print("="*88)
print(f"{'window':>14} {'L':>4} {'#in-zeros':>10} {'N':>4} {'cond':>6} {'margin':>13} {'next_out_zero':>13}")
for L in [20,30,40,50,60,70]:
    a=14.0; b=a+L; W=1.0; nz=int(((ZT>=a)&(ZT<=b)).sum()); N=max(4,round(0.75*L*2*W/np.pi))
    outz=ZT[ZT>b]; nxt=float(outz.min()) if len(outz) else np.nan
    mg,cond=margin(a,b,W,N)
    print(f"{f'[{a:.0f},{b:.0f}]':>14} {L:>4} {nz:>10} {N:>4} {cond:>6.1f} {mg:>+13.4e} {nxt:>13.3f}",flush=True)
print("="*88)
print("If margin thins as L (and #in-zeros) grows -> confirms §858: more captured zeros = residual from further tails.")
print("="*88)

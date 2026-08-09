#!/usr/bin/env python3
"""§892: height-ladder of the scalar-criterion quantities. Fixed 4W=4, frac=0.65, windows climbing height 15→1000.
Track μ (depth of the single negative eigenvalue of M_σ=arch−prime — the PRIME-defect/RH quantity), p (pole
strength, explicit/archimedean), (v·w)² (alignment), π_w=p(v·w)², and margin. Where does μ cross 0 (corner→robust)?
Does μ obey a clean law vs floor/height? The pole side (p) is explicit; μ carries the RH content — its law is the goal."""
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
NMAX=5000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)

def run(a,b,fourW,frac=0.65,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RGRID=np.arange(-(b+60),(b+60),dr)
    OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    Msw=B@(arch-prime)@B; polew=B@pole@B
    ev,V=eigh(Msw); mu=-float(ev[0]); w=V[:,0]; nneg=int((ev<0).sum())
    pev,pvec=eigh(polew); p=float(pev[-1]); v=pvec[:,-1]
    floor=float(eigh(B@arch@B+polew,eigvals_only=True)[0])
    margin=float(eigh(Msw+polew,eigvals_only=True)[0])
    return N,nneg,mu,p,float(v@w)**2,p*float(v@w)**2,floor,margin

print("="*120)
print("§892 height-ladder (4W=4): μ=defect depth (RH quantity), p=pole strength (explicit), (v·w)², π_w, floor, margin")
print("="*120)
print(f"{'height':>8} {'N':>3} {'#neg':>4} {'μ (defect depth)':>16} {'p (pole)':>10} {'(v·w)²':>8} {'π_w=p(v·w)²':>12} {'floor':>7} {'margin':>12}")
mids=[]
for (a,b) in [(14,20),(17,26),(22,34),(30,45),(40,60),(55,80),(75,110),(110,150),(160,210),(230,300),(400,500),(700,820)]:
    N,nneg,mu,p,vw2,pi_w,floor,margin=run(a,b,4.0)
    h=(a+b)/2; mids.append((h,mu,p,vw2,pi_w,floor))
    print(f"{h:>8.0f} {N:>3} {nneg:>4} {mu:>+16.4e} {p:>10.3e} {vw2:>8.4f} {pi_w:>12.3e} {floor:>7.4f} {margin:>+12.3e}")
print("="*120)
# fit μ (where positive=corner) vs floor and vs height
arr=np.array(mids); corner=arr[arr[:,1]>0]
if len(corner)>=3:
    h,mu,p,vw2,pi_w,fl=corner.T
    b1,b0=np.polyfit(np.log(h),np.log(mu),1); print(f"corner μ ~ height^{b1:.2f} (log-log slope), intercept {b0:.2f}")
    b1f,b0f=np.polyfit(fl,np.log(mu),1); print(f"corner log μ ~ {b1f:.2f}·floor + {b0f:.2f}")
print("μ crosses 0 at the corner→robust transition. p (pole) is explicit/archimedean; μ carries the RH content.")
print("If μ has a clean law (vs floor/height/prime-gap) → the arithmetic origin of the razor.")
print("="*120)

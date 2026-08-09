#!/usr/bin/env python3
"""§912: is the RAZOR (π_w/μ ≈ 1) bottom-bounded? For each 4W, walk corner windows from the bottom upward and record
π_w/μ. If ratio ≈ 1 only near the bottom (h<~40) and GROWS with height (pole over-covers upward), the genuinely-hard
razor is confined to a bounded bottom region — the refined finite-proof path (§911)."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh, pinvh
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3200; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def ratio(a,b,fourW,frac=0.85,dr=0.005):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+50),(b+50),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    gw=np.linalg.eigvalsh(Gm)
    if gw[0]<=1e-13: return None
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw2,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw2))@gv.T
    Msw=B@(arch-prime)@B; polew=B@pole@B; ev,V=eigh(Msw)
    if ev[0]>=0: return (N,None,None)   # no corner (PSD)
    mu=-ev[0]; w=V[:,0]; pev,pvec=eigh(polew); p=pev[-1]; v=pvec[:,-1]
    return (N,mu,p*(v@w)**2/mu)
print("="*96)
print("§912 is the razor bottom-bounded? π_w/μ across corner windows per 4W (climbing height)")
print("="*96)
WINS=[(14,30),(14,44),(24,44),(30,50),(30,60),(40,64),(45,70),(55,80),(60,90),(75,105)]
print(f"{'4W':>4} | "+" ".join(f"h{(a+b)//2}" for a,b in WINS))
for fourW in [4.0,5.0,6.0,7.0,8.0]:
    row=[]
    for (a,b) in WINS:
        r=ratio(a,b,fourW)
        if r is None: row.append("  cnd "); continue
        N,mu,rt=r
        if mu is None: row.append("  PSD "); continue
        row.append(f"{rt:5.2f}")
    print(f"{fourW:>4.1f} | "+" ".join(row))
print("="*96)
print("PSD = no corner (M_σ⪰0). ratio π_w/μ: ~1.0x = RAZOR (hard); >1.5 = pole over-covers (easy).")
print("If razor (~1.0) only at low h and ratio GROWS with height → hard part bottom-bounded → refined finite path.")
print("="*96)

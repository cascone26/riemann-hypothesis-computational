#!/usr/bin/env python3
"""§890: is #neg(M_σ)=#negative eigenvalues of the exact symbol form arch−prime always ≤ rank(pole)=2, and does
the pole always lift them to ≥0? Sweep deeper corners (4W up to 8), many windows/heights, near & below G-capacity.
If #neg(M_σ)≤2 everywhere → clean finite necessary structure (pole is a rank-2 lifter). If it grows → band-limited
coupling (Schur) is essential beyond the pole. Report #neg, the pole's 2 eigenvalues, and full-margin sign."""
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
NMAX=8000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)

def run(a,b,fourW,frac,dr=0.01):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RGRID=np.arange(-(b+40),(b+40),dr)
    OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr
    gw=np.linalg.eigvalsh(G)
    if gw[0]<=1e-13: return N,Ncap,None,None,None,None  # G not PD, skip
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gv=np.linalg.eigh(G)[1]; B=gv@np.diag(1/np.sqrt(gw))@gv.T
    Msw=B@(arch-prime)@B; polew=B@pole@B
    ev=eigh(Msw,eigvals_only=True); nneg=int((ev<-1e-12).sum())
    pev=eigh(polew,eigvals_only=True)[::-1]
    mfull=float(eigh(Msw+polew,eigvals_only=True)[0])
    return N,Ncap,nneg,float(ev[0]),(pev[0],pev[1]),mfull

print("="*114)
print("§890 #neg eigenvalues of exact symbol form M_σ=arch−prime, vs pole rank(=2). Does pole always lift to ≥0?")
print("="*114)
print(f"{'window':>12} {'4W':>4} {'frac':>4} {'N':>3} {'Ncap':>5} {'#neg(M_σ)':>9} {'min M_σ':>11} {'pole eig1,eig2':>20} {'full margin':>12}")
maxneg=0; anyfail=False
cfgs=[]
for (a,b) in [(14,44),(14,34),(30,60),(60,90),(120,150)]:
    for fourW in [4.0,5.0,6.0,7.0,8.0]:
        for frac in [0.7,0.85,0.95]:
            cfgs.append((a,b,fourW,frac))
for (a,b,fourW,frac) in cfgs:
    N,Ncap,nneg,mmin,pev,mfull=run(a,b,fourW,frac)
    if nneg is None: continue
    maxneg=max(maxneg,nneg); flag=""
    if nneg>2: flag=" <-- #neg>2!"; anyfail=True
    if mfull<-1e-11: flag+=" NEGmargin"
    if nneg>0:  # only print the interesting (corner) rows + a few robust
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {frac:>4.2f} {N:>3} {Ncap:>5.1f} {nneg:>9} {mmin:>+11.3e} {f'{pev[0]:.2e},{pev[1]:.1e}':>20} {mfull:>+12.3e}{flag}")
print("="*114)
print(f"MAX #neg(M_σ) over all {len(cfgs)} configs = {maxneg}.  Any #neg>2 or negative full margin: {anyfail}")
print("If max #neg ≤ 2 everywhere and full margin ≥0 → rank-2 pole always suffices to lift M_σ (clean necessary structure).")
print("="*114)

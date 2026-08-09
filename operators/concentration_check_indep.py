#!/usr/bin/env python3
"""§888 (independent): concentration operator K = B_W·1_{σ<0}·B_W. λ_max(K)=max fraction of |φ|² mass a
bandwidth-4W φ can put on {σ<0}. Naive Szegő bound: margin ≳ floor − λ_max·|min σ|? Track across heights.
Independent of the workflow's operators/concentration_operator.py (cross-check by independence)."""
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

def run(a,b,fourW,frac=0.7,dr=0.002):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(frac*Ncap))
    RGRID=np.arange(-(b+120),(b+120),dr)
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
    margin=float(eigh(arch+pole-prime,G,eigvals_only=True)[0])
    # symbol on grid
    Pw=(1/np.pi)*(lw[:,None]*np.cos(np.outer(ln,RGRID))).sum(axis=0); sig=OMEGA/(2*np.pi)-Pw
    negmask=sig<0.0
    floor=float(eigh(arch+pole,G,eigvals_only=True)[0])
    minsig=float(sig[(RGRID>=a-2)&(RGRID<=b+2)].min())
    # concentration operator: M_ij = ∫ P_i 1_{σ<0} P_j ; K=G^{-1}M ; λ_max via gen-eig(M,G)
    Pn=P*negmask[None,:]; M=(Pn@P.T)*dr; M=0.5*(M+M.T)
    lam=eigh(M,G,eigvals_only=True); lammax=float(lam[-1])
    # also: the actual min of ∫|φ|²σ / ‖φ‖² (symbol/Szegő margin) via gen-eig(Msig, G)
    Ms=(P*sig[None,:])@P.T*dr; Ms=0.5*(Ms+Ms.T); szmargin=float(eigh(Ms,G,eigvals_only=True)[0])
    return N,Ncap,margin,floor,minsig,lammax,szmargin

print("="*118)
print("§888 concentration operator λ_max(K) and the Szegő-symbol margin vs the true Weil margin, across heights (4W=4)")
print("="*118)
print(f"{'window':>12} {'N':>3} {'Ncap':>5} {'true margin':>11} {'floor':>7} {'min σ':>8} {'λmax(K)':>8} {'floor−λmax|minσ|':>16} {'Szegő margin':>12}")
for (a,b) in [(14,44),(60,90),(120,150),(300,330),(1000,1030)]:
    N,Ncap,margin,floor,minsig,lammax,szm=run(a,b,4.0)
    naive=floor-lammax*abs(minsig)
    print(f"{f'[{a},{b}]':>12} {N:>3} {Ncap:>5.1f} {margin:>+11.3e} {floor:>7.4f} {minsig:>+8.3f} {lammax:>8.4f} {naive:>+16.4f} {szm:>+12.4e}")
print("="*118)
print("λmax(K)=max negative-set mass a 4W-band φ can capture. If Szegő margin ≈ true margin → symbol frame VALID.")
print("naive bound floor−λmax|minσ| is the crude Szegő lower bound (likely negative=too weak; the real content is finer).")
print("="*118)

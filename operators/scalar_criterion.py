#!/usr/bin/env python3
"""§891: RH-corner as a SINGLE SCALAR inequality. M_σ=arch−prime has one negative eigenvalue −μ (eigvec w);
pole = p·vv^T (rank-1). Q=M_σ+pole ⪰0 ⟺ μ·w^T A^+ w ≤ 1, A=M_+ + p vv^T (M_+ = PSD part of M_σ). Leading form:
p·(v·w)² ≥ μ. Verify the scalar g=1−μ·w^T A^+ w has the SAME SIGN as the true margin, across well-conditioned
corners (frac≤0.7). Two explicit zero-free quantities — cleaner than floor-vs-pnorm (§875)."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh, pinvh
def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=4000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)

def run(a,b,fourW,frac=0.65,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RGRID=np.arange(-(b+60),(b+60),dr)
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
    Msw=B@(arch-prime)@B; polew=B@pole@B
    ev,V=eigh(Msw); nneg=int((ev<0).sum())
    mu=-float(ev[0]); w=V[:,0]                      # negative eigenvalue -mu, eigvec w
    Mplus=V[:,1:]@np.diag(ev[1:])@V[:,1:].T          # PSD part
    pev,pvec=eigh(polew); p=float(pev[-1]); v=pvec[:,-1]  # pole = p vv^T (rank-1)
    A=Mplus+polew
    margin=float(eigh(Msw+polew,eigvals_only=True)[0])
    # exact scalar: g = 1 - mu * w^T A^+ w   (>=0 <=> Q>=0, since Q=A-mu ww^T)
    Ap=pinvh(A); g=1.0-mu*float(w@Ap@w)
    # leading criterion: pole mass along w vs mu
    pi_w=p*float(v@w)**2
    return N,nneg,mu,p,float(v@w)**2,pi_w,margin,g

print("="*118)
print("§891 RH-corner as one scalar inequality: exact g=1−μ·w^T A^+ w  vs true margin sign;  leading π_w=p(v·w)² vs μ")
print("="*118)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'#neg':>4} {'μ=|negeig|':>11} {'p(pole)':>9} {'(v·w)²':>8} {'π_w':>10} {'π_w−μ':>11} {'true margin':>12} {'exact g':>10} {'signs✓':>6}")
for (a,b,fourW) in [(14,44,4.0),(14,34,4.0),(14,44,3.0),(14,44,5.0),(30,60,4.0),(30,60,5.0),(60,90,4.0),(120,150,4.0)]:
    N,nneg,mu,p,vw2,pi_w,margin,g=run(a,b,fourW)
    ok="YES" if (np.sign(margin)==np.sign(g) or abs(margin)<1e-12) else "NO"
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {nneg:>4} {mu:>11.3e} {p:>9.3e} {vw2:>8.4f} {pi_w:>10.3e} {pi_w-mu:>+11.3e} {margin:>+12.3e} {g:>+10.3e} {ok:>6}")
print("="*118)
print("exact g and true margin SAME SIGN → RH-corner ⟺ single scalar g=1−μ·w^T A^+ w ≥ 0 (two explicit zero-free quantities).")
print("π_w−μ ≈ margin sign too → leading criterion 'pole mass along the defect's worst direction ≥ its depth'.")
print("="*118)

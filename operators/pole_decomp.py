#!/usr/bin/env python3
"""§889: CLEAN DECOMPOSITION. Is arch−prime EXACTLY the symbol compression M_σ (σ=Ω/2π−P)? Then Q = M_σ + pole,
pole = rank-≤2 PSD outer product. Test at the corner: (1) ‖(arch−prime) − M_σ‖ (should be ~0 → exact symbol),
(2) margin(arch−prime) [no pole] vs margin(arch+pole−prime) [full] — does the rank-2 pole flip the sign?,
(3) #negative eigenvalues of M_σ vs G, and how well the pole's range covers that negative eigenspace.
Resolves §885(pole negligible) vs §888(pole essential): pole is negligible in the COUPLING but may be essential on
the DIAGONAL/margin. HONEST arithmetic, no interpretation until the numbers are in."""
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
    # symbol compression M_sig = compression of mult by (Omega/2pi - P(r))
    Pw=(1/np.pi)*(lw[:,None]*np.cos(np.outer(ln,RGRID))).sum(axis=0); sig=OMEGA/(2*np.pi)-Pw
    Msig=(P*sig[None,:])@P.T*dr; Msig=0.5*(Msig+Msig.T)
    # (1) is arch-prime == Msig ?
    AmP=arch-prime; sym_err=np.linalg.norm(AmP-Msig)/np.linalg.norm(AmP)
    # whiten to compare eigen-spectra vs G
    gw,gv=eigh(G); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    AmPw=B@AmP@B; polew=B@pole@B; Qw=AmPw+polew
    ev_AmP=eigh(AmPw,eigvals_only=True); ev_Q=eigh(Qw,eigvals_only=True)
    m_nopole=float(ev_AmP[0]); m_full=float(ev_Q[0])
    nneg=int((ev_AmP<0).sum()); poler=int((eigh(polew,eigvals_only=True)>1e-10).sum())
    # (3) how well does pole's range cover the negative eigenspace of AmP?
    negvecs=eigh(AmPw)[1][:,ev_AmP<0]  # columns = negative eigvecs (whitened)
    pw_ev,pw_vec=eigh(polew); prange=pw_vec[:,pw_ev>1e-10]  # pole range
    # projection of each neg eigenvector onto pole range
    cov=np.linalg.norm(prange.T@negvecs,axis=0) if negvecs.shape[1]>0 and prange.shape[1]>0 else np.array([])
    return N,Ncap,sym_err,m_nopole,m_full,nneg,poler,ev_AmP[:min(4,N)],cov

print("="*112)
print("§889 Q = M_σ (exact symbol, arch−prime) + pole (rank≤2 PSD). Does the pole flip the corner margin sign?")
print("="*112)
for (a,b,fourW) in [(14,44,4.0),(60,90,4.0),(120,150,4.0)]:
    N,Ncap,se,mnp,mf,nneg,pr,evs,cov=run(a,b,fourW)
    print(f"\n[{a},{b}] 4W={fourW} N={N} (Ncap={Ncap:.1f}):")
    print(f"  ‖(arch−prime) − M_σ‖/‖·‖ = {se:.2e}   → arch−prime {'IS' if se<1e-6 else 'is NOT'} the exact symbol compression")
    print(f"  margin WITHOUT pole (M_σ) = {mnp:+.4e}    margin WITH pole (full Q) = {mf:+.4e}    pole rank={pr}")
    print(f"  #negative eigenvalues of M_σ = {nneg} ;  lowest M_σ eigenvalues: {np.array2string(evs,precision=4,suppress_small=False)}")
    print(f"  pole-range coverage of M_σ negative eigvec(s) (‖proj‖, 1=fully covered): {np.array2string(cov,precision=3)}")
print("="*112)
print("If margin flips −→+ with the rank≤2 pole AND pole-range covers M_σ's negative eigenspace → RH-corner ⟺")
print("'a rank-2 PSD pole term lifts the indefinite exact-symbol form M_σ to PSD'. Clean finite reformulation.")
print("="*112)

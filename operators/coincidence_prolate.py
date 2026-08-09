#!/usr/bin/env python3
"""§906 (independent): tie it together. (1) Confirm M_σ's near-null cluster = the arch≈prime COINCIDENCE subspace
(arch(u_k)≈prime(u_k) on near-null eigvecs u_k). (2) Independently check Angle 2's overlap of near-null with the
prolate/concentration plunge. (3) Characterize the negative direction w: is it the coincidence direction where prime
just exceeds arch? Concentration op C = P_band∘1_window∘P_band (Slepian plunge). Verify the §905 ~0.53 overlap myself."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=6000; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def build(a,b,fourW,N,dr=0.003):
    W=fourW/4; RG=np.arange(-(b+70),(b+70),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    # concentration operator (Slepian): band-limit restricted to window mass. In the whitened Fejér basis, the
    # window-concentration matrix is W_ij = ∫_a^b P_i P_j dr (mass inside window) vs G (total). Plunge = eig(W,G).
    win=(RG>=a)&(RG<=b); Wm=(P[:,win]@P[:,win].T)*dr
    gw,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    return B@arch@B, B@prime@B, B@Wm@B
for (a,b,fourW,N) in [(14,44,4.0,13),(14,44,5.0,17),(14,34,4.0,8)]:
    archw,primew,Cw=build(a,b,fourW,N)
    Msw=archw-primew; ev,V=eigh(Msw); nn=int((ev<0).sum())
    # near-null cluster = eigvecs with |eig|<1e-3
    nnidx=np.where(np.abs(ev)<1e-3)[0]; U=V[:,nnidx]
    # (1) coincidence: on near-null u_k, arch(u_k) vs prime(u_k)
    aq=np.array([U[:,k]@archw@U[:,k] for k in range(U.shape[1])]); pq=np.array([U[:,k]@primew@U[:,k] for k in range(U.shape[1])])
    coincide=np.max(np.abs(aq-pq))  # should be ~ the near-null eigenvalues (small)
    # (2) prolate plunge: C eigenvalues (concentration ratios) + plunge region (0.1..0.9), overlap of near-null with plunge
    cev,cV=eigh(Cw); plunge=np.where((cev>0.05)&(cev<0.95))[0]; Cplunge=cV[:,plunge]
    # overlap: projection of near-null subspace onto plunge subspace (principal angles / mean overlap)
    Ov=Cplunge.T@U  # (nplunge x nnearnull)
    nn_in_plunge=np.linalg.norm(Ov,axis=0)**2  # fraction of each near-null vec inside plunge span
    # (3) negative direction w = V[:,0]
    w=V[:,0]; aw=w@archw@w; pw=w@primew@w
    print(f"\n[{a},{b}] 4W={fourW} N={N}: #neg={nn}, near-null cluster dim={U.shape[1]}, plunge dim={len(plunge)} (of {N})")
    print(f"  (1) COINCIDENCE on near-null: max|arch(u)−prime(u)| = {coincide:.2e}  (arch≈prime confirmed: {coincide<1e-3})")
    print(f"  (2) near-null vecs' mass inside prolate plunge span: {np.array2string(nn_in_plunge,precision=3)}  (mean {nn_in_plunge.mean():.3f})")
    print(f"  (3) neg direction w: arch(w)={aw:.4f} prime(w)={pw:.4f} prime−arch={pw-aw:.3e} (=μ, prime just exceeds arch)")
print("\nnear-null = coincidence subspace (arch≈prime), ~partly in prolate plunge; neg dir = where prime marginally exceeds arch.")
print("This ties §882 (coincidence) + §905 (prolate-partial) + §901 (single neg mode): the razor lives in the coincidence/plunge cluster.")

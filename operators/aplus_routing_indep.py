#!/usr/bin/env python3
"""§901 (independent): decompose the razor quantity w^T A^+ w = Σ_k (w·a_k)²/α_k over A-eigenmodes (A=M_+ + p vv^T).
Which modes drive it — the near-null (small α_k), the pole (α_k≈p), or bulk? Pins down exactly which modes make the
exact criterion g=1−μ w^T A^+ w razor. Independent of the workflow's operators/aplus_routing.py (cross-check)."""
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
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    return B@(arch-prime)@B, B@pole@B
print("="*100)
print("§901 A^+ routing: which A-eigenmodes drive w^T A^+ w (the razor)? A=M_+ + p vv^T, M_+=PSD part of M_σ")
print("="*100)
for (a,b,fourW,N) in [(14,44,4.0,13),(14,44,5.0,17),(14,34,4.0,8)]:
    Msw,polew=build(a,b,fourW,N)
    ev,V=eigh(Msw); mu=-ev[0]; w=V[:,0]
    Mplus=V[:,1:]@np.diag(ev[1:])@V[:,1:].T
    A=Mplus+polew
    al,a_k=eigh(A)  # A eigenpairs ascending
    contrib=np.array([ (w@a_k[:,k])**2/al[k] for k in range(N) ])
    wAw=contrib.sum(); g=1-mu*wAw
    margin=float(eigh(Msw+polew,eigvals_only=True)[0])
    p=float(eigh(polew,eigvals_only=True)[-1])
    # rank contributions
    order=np.argsort(contrib)[::-1]
    print(f"\n[{a},{b}] 4W={fourW} N={N}: μ={mu:.3e} p(pole)={p:.3e} w^T A^+ w={wAw:.3e} g={g:+.3e} margin={margin:+.3e}")
    print(f"  A eigenvalues (asc): {np.array2string(al[:4],precision=2)} ... {np.array2string(al[-2:],precision=2)}")
    print(f"  top-4 razor contributions (w·a_k)²/α_k: "+", ".join(f"α={al[order[i]]:.2e}→{contrib[order[i]]:.2e}({100*contrib[order[i]]/wAw:.0f}%)" for i in range(4)))
    # how much from the smallest-alpha (near-null) half vs pole (largest alpha)?
    near=contrib[al<1e-3].sum(); polec=contrib[al>0.5*p].sum() if p>1e-8 else 0
    print(f"  from near-null A-modes (α<1e-3): {100*near/wAw:.1f}%   from pole-scale modes (α>p/2): {100*polec/wAw:.1f}%")
print("="*100)
print("If w^T A^+ w is dominated by SMALL-α (near-null) A-modes → the razor IS the near-null coupling (§900 confirmed).")
print("="*100)

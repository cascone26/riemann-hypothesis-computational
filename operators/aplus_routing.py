#!/usr/bin/env python3
"""ANGLE B — How does A^+ route the razor? Decompose w^T A^+ w = Σ_k (w·a_k)² / α_k
and identify which modes (near-null, pole, bulk) dominate the RH-corner razor."""
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

def run_angle_b(a,b,fourW,frac=0.65,dr=0.004):
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

    # Eigendecompose A = Σ_k α_k a_k a_k^T
    alpha_k, a_k = eigh(A)

    # Compute w·a_k for each mode
    w_dot_ak = w @ a_k  # shape (N,)

    # Compute per-mode contributions: (w·a_k)² / α_k
    contrib = np.zeros(N)
    for i in range(N):
        if abs(alpha_k[i]) > 1e-14:
            contrib[i] = (w_dot_ak[i]**2) / alpha_k[i]
        else:
            contrib[i] = 0

    # Total w^T A^+ w
    wAw = np.sum(contrib)

    # Compute g = 1 - μ·w^T A^+ w
    g = 1.0 - mu * wAw

    # Classify modes: near-null, pole, bulk
    pole_idx = np.argmax(alpha_k)
    null_thresh = 1e-5  # modes with α_k < 1e-5 are "near-null"
    null_mask = alpha_k < null_thresh
    bulk_mask = (alpha_k >= null_thresh) & (alpha_k < alpha_k[pole_idx] * 0.9)

    # Sort modes by contribution descending
    sort_idx = np.argsort(-contrib)

    # Sum contributions by class
    contrib_null = np.sum(contrib[null_mask]) if np.any(null_mask) else 0.0
    contrib_pole = contrib[pole_idx]
    contrib_bulk = np.sum(contrib[bulk_mask]) if np.any(bulk_mask) else 0.0

    return {
        'N': N, 'a': a, 'b': b, 'fourW': fourW,
        'mu': mu, 'p': p, 'margin': margin, 'g': g,
        'alpha_k': alpha_k, 'w_dot_ak': w_dot_ak, 'contrib': contrib,
        'wAw': wAw, 'sort_idx': sort_idx, 'pole_idx': pole_idx,
        'null_mask': null_mask, 'bulk_mask': bulk_mask,
        'contrib_null': contrib_null, 'contrib_pole': contrib_pole, 'contrib_bulk': contrib_bulk,
        'margin_sign': np.sign(margin), 'g_sign': np.sign(g),
    }

if __name__ == '__main__':
    print("="*120)
    print("ANGLE B: A^+ routing analysis — which modes dominate w^T A^+ w = Σ_k (w·a_k)² / α_k?")
    print("="*120)

    corners = [(14,44,4.0), (14,44,5.0)]

    for a,b,fourW in corners:
        result = run_angle_b(a,b,fourW)
        print(f"\n[{a},{b}] 4W={fourW} N={result['N']}")
        print(f"  μ={result['mu']:.3e}, p={result['p']:.3e}, margin={result['margin']:+.3e}, g={result['g']:+.3e}")
        print(f"  margin_sign={result['margin_sign']:+.0f}, g_sign={result['g_sign']:+.0f} → signs match: {result['margin_sign']==result['g_sign']}")
        print(f"\n  w^T A^+ w decomposition by mode (top 10 contributions):")
        print(f"  {'rank':>4} {'i':>4} {'α_k':>12} {'w·a_k':>12} {'contrib':>12} {'class':>8}")
        for rank, i in enumerate(result['sort_idx'][:10]):
            if result['null_mask'][i]:
                mode_class = 'null'
            elif i == result['pole_idx']:
                mode_class = 'POLE'
            elif result['bulk_mask'][i]:
                mode_class = 'bulk'
            else:
                mode_class = '?'
            print(f"  {rank+1:>4} {i:>4} {result['alpha_k'][i]:>12.3e} {result['w_dot_ak'][i]:>12.3e} {result['contrib'][i]:>12.3e} {mode_class:>8}")

        total = result['wAw']
        print(f"\n  Total w^T A^+ w = {total:.6e}")
        print(f"  Contributions by class:")
        print(f"    Near-null modes: {result['contrib_null']:.6e} ({result['contrib_null']/total*100 if total>0 else 0:.2f}%)")
        print(f"    Pole mode:       {result['contrib_pole']:.6e} ({result['contrib_pole']/total*100 if total>0 else 0:.2f}%)")
        print(f"    Bulk modes:      {result['contrib_bulk']:.6e} ({result['contrib_bulk']/total*100 if total>0 else 0:.2f}%)")

    print("="*120)

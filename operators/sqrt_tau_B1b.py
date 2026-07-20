#!/usr/bin/env python3
"""B1b: the CORRECT SUSY test — is the SYMMETRIC square root sqrt(tau) local (quasi-banded), c-stably?

B1 tested the triangular Cholesky factor (delocalized). But the SUSY lead A12/A13 was about the SYMMETRIC
root tau^(1/2) (98% band energy), a different factorization. tau = R R with R = tau^(1/2) = V diag(sqrt L) V^T.
sqrt(tau) is dominated by the O(1) eigenvalues, so double-precision eigh of the Arb-assembled tau suffices
for a locality measure (the tiny eigenvalue contributes ~1e-11, negligible). If sqrt(tau) is banded and
c-stable, Vector B (manifest local factorization) stays LIVE.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import numpy as np
from flint import arb, ctx
from arb_ldlt_certify import build_arb_tau

def tau_np(c,N,prec=2000):
    ctx.prec=prec; A,DIM=build_arb_tau(c,N,prec)
    T=np.zeros((DIM,DIM))
    for i in range(DIM):
        for j in range(i,DIM):
            v=float(A[i,j].mid()); T[i,j]=v; T[j,i]=v
    return T,DIM

def band_energy(M):
    DIM=M.shape[0]; tot=float(np.sum(M*M)) or 1.0
    ii,jj=np.indices((DIM,DIM)); b=np.abs(ii-jj)
    return {k: float(np.sum((M*M)[b<=k]))/tot for k in (1,2,3,5,8,12,20)}

def main():
    N=64
    print(f"=== B1b: locality of symmetric sqrt(tau) (N={N}) ===")
    out=[]
    for c in [14,50,200,800]:
        T,DIM=tau_np(c,N)
        w,V=np.linalg.eigh(T)
        w=np.clip(w,0,None)
        R=(V*np.sqrt(w)[None,:])@V.T          # symmetric sqrt
        be=band_energy(R)
        prof=[float(np.mean(np.abs(np.diagonal(R,-r)))) for r in range(0,10)]
        print(f"\nc={c:4d}: sqrt(tau) band energy within |i-j|<=k:")
        print("   " + "  ".join(f"k{k}:{100*be[k]:.0f}%" for k in (1,2,3,5,8,12,20)))
        print(f"   subdiag |R| r=0..6: " + " ".join(f"{prof[r]:.2e}" for r in range(7)))
        out.append(dict(c=c, band_energy=be, subdiag=prof))
    json.dump(out,open("results/sqrt_tau_B1b.json","w"),indent=2)
    print("\n=== VERDICT (B1b) ===")
    k3=[r['band_energy'][3] for r in out]
    banded=all(x>0.75 for x in k3); stable=(max(k3)-min(k3))<0.15
    if banded and stable:
        print(f"  => sqrt(tau) is LOCAL & c-stable ({[f'{100*x:.0f}%' for x in k3]}). Vector B LIVE:")
        print("     a local symmetric factorization root exists — pursue its c->inf limit form.")
    elif banded:
        print(f"  => sqrt(tau) banded but drifts ({[f'{100*x:.0f}%' for x in k3]}); check limit.")
    else:
        print(f"  => sqrt(tau) not local ({[f'{100*x:.0f}%' for x in k3]}); SUSY-locality lead does not persist.")
    print("wrote results/sqrt_tau_B1b.json")

if __name__=="__main__":
    main()

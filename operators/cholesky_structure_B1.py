#!/usr/bin/env python3
"""VECTOR B1: structure of the Cholesky/LDL^T factor of tau. Does tau = M M^T with M LOCAL (banded)
and structurally STABLE across the prime cutoff c? (SUSY lead A12/A13: tau^1/2 is quasi-local, 98% band
energy.) If M is banded with a c-stable profile, a manifest-positivity factorization limit is plausible
(Vector B); if M is dense/unstructured, factorization is non-informative and B is parked.

Computes L (unit lower) and pivots D via Arb LDL^T, forms M = L * diag(sqrt(D)) (real since tau PD), and
measures band energy concentration and off-diagonal decay, across c.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import numpy as np
from flint import arb, ctx
from arb_ldlt_certify import build_arb_tau

def ldlt_np(c, N, prec=2500):
    ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec)
    d=[None]*DIM; Lf=[[arb(0)]*DIM for _ in range(DIM)]
    for i in range(DIM):
        s=A[i,i]
        for k in range(i): s=s-Lf[i][k]*Lf[i][k]*d[k]
        d[i]=s
        for j in range(i+1,DIM):
            t=A[j,i]
            for k in range(i): t=t-Lf[j][k]*Lf[i][k]*d[k]
            Lf[j][i]=t/d[i]
    L=np.eye(DIM)
    for j in range(DIM):
        for i in range(j+1,DIM):
            L[i,j]=float(Lf[i][j].mid())
    dd=np.array([max(float(x.mid()),0.0) for x in d])
    M=L*np.sqrt(dd)[None,:]   # M = L diag(sqrt d); tau = M M^T
    return M, DIM

def band_energy(M):
    """fraction of Frobenius^2 energy within |i-j|<=k, for several k."""
    DIM=M.shape[0]; tot=float(np.sum(M*M)) or 1.0
    ii,jj=np.indices((DIM,DIM)); band=np.abs(ii-jj)
    return {k: float(np.sum((M*M)[band<=k]))/tot for k in (0,1,2,3,5,8,12,20)}

def main():
    N=64
    print(f"=== B1: Cholesky factor locality of tau=M M^T (N={N}) ===")
    out=[]
    for c in [14,50,200]:
        M,DIM=ldlt_np(c,N)
        be=band_energy(M)
        # off-diagonal decay: mean |M[i,i-r]| vs r (subdiagonal profile)
        prof=[float(np.mean(np.abs(np.diagonal(M,-r)))) for r in range(0,15)]
        print(f"\nc={c:4d}: band-energy within |i-j|<=k:")
        print("   " + "  ".join(f"k{k}:{100*be[k]:.0f}%" for k in (0,1,2,3,5,8,12,20)))
        print(f"   subdiagonal |M| profile r=0..8: " + " ".join(f"{prof[r]:.2e}" for r in range(9)))
        out.append(dict(c=c, band_energy=be, subdiag_profile=prof))
    json.dump(out, open("results/cholesky_structure_B1.json","w"), indent=2)
    # verdict
    print("\n=== VERDICT (B1) ===")
    k5=[r['band_energy']['5'] if '5' in r['band_energy'] else r['band_energy'][5] for r in out]
    # json keys become strings on reload but here in-memory keys are ints
    k5=[r['band_energy'][5] for r in out]
    banded = all(x>0.8 for x in k5)
    stable = (max(k5)-min(k5))<0.15
    if banded and stable:
        print(f"  => M is LOCAL (>80% energy within 5 bands) and c-STABLE ({[f'{100*x:.0f}%' for x in k5]}).")
        print("     Manifest-positivity factorization limit is plausible — Vector B stays LIVE.")
    elif banded:
        print(f"  => M banded but c-profile drifts ({[f'{100*x:.0f}%' for x in k5]}); investigate limit.")
    else:
        print(f"  => M is NOT local ({[f'{100*x:.0f}%' for x in k5]}); factorization non-informative, park B.")
    print("wrote results/cholesky_structure_B1.json")

if __name__=="__main__":
    main()

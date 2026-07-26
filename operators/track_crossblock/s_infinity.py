#!/usr/bin/env python3
"""P2 — the construction question on S∞. τ⪰0 ⟺ low-freq Schur complement S ⪰ 0 (A9). S converges to a fixed
finite S∞ as the high-freq bandwidth →∞. The reset target: does S∞ carry a finite-dimensional POLARIZATION
(a structure making its positivity manifest — Gram/local/arithmetic), the way function-field H¹ does? Or is it
degenerate/nonlocal like every other view?

Tests, for Groskin τ, low-freq block L={|n|≤M}, high-freq H={|n|>M} Schur-eliminated, N (=high bandwidth) growing:
 (1) CONVERGENCE: does S(N) → S∞ as N grows? (A9 claim)
 (2) KERNEL + GAP: dim of near-kernel of S∞, and the GAP above it (smallest nonzero eig) — the quotient-margin.
 (3) UNIFORMITY in c: does gap(S∞) stay bounded below, or →0 as c grows? (uniform => min-max attackable)
 (4) POLARIZATION/LOCALITY: is S∞ (or its non-kernel part) recognizable — banded/Toeplitz/low-rank-arithmetic,
     i.e. does it have a manifest-positive structure, or is it dense/nonlocal (= the wall)?
Float64 for structure (the gap ~O(1e-2..1)? or tiny? — we'll see; the O(1) structure is reliable).
"""
import importlib.util, os, numpy as np
GPATH=os.path.expanduser('~/projects/riemann/external/groskin_2607.02828/code/arb_ldlt_certify.py')
spec=importlib.util.spec_from_file_location('gc',GPATH); g=importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
def tof(x):
    try: return float(x)
    except Exception:
        s=str(x)
        if '+/-' in s: s=s.split('+/-')[0].strip('[] ')
        return float(s)
def schur(c,N,M,prec=1000):
    A,DIM=g.build_arb_tau(c,N,prec)
    T=np.array([[tof(A[i,j]) for j in range(DIM)] for i in range(DIM)])
    idx=np.arange(DIM); n=idx-N
    Lm=np.where(np.abs(n)<=M)[0]; Hm=np.where(np.abs(n)>M)[0]
    TLL=T[np.ix_(Lm,Lm)]; TLH=T[np.ix_(Lm,Hm)]; THH=T[np.ix_(Hm,Hm)]
    S=TLL - TLH@np.linalg.solve(THH,TLH.T)
    return S

M=10
print(f"P2 S∞ construction test. Low-freq block |n|<= {M} (dim {2*M+1}). High-freq Schur-eliminated, N growing.")
for c in [13,29,50]:
    print(f"\n===== c={c} =====")
    Ss={}
    for N in [15,25,35]:
        Ss[N]=schur(c,N,M)
    # (1) convergence
    d1=np.linalg.norm(Ss[25]-Ss[15])/np.linalg.norm(Ss[15]); d2=np.linalg.norm(Ss[35]-Ss[25])/np.linalg.norm(Ss[25])
    print(f" (1) convergence: ||S25-S15||/||S15||={d1:.3e}, ||S35-S25||/||S25||={d2:.3e}  (→0 = S∞ converges)")
    Sinf=Ss[35]
    w=np.sort(np.linalg.eigvalsh(Sinf))
    # (2) kernel + gap: count near-zero eigs, gap above
    tol=1e-6*np.max(np.abs(w)); ker=int(np.sum(np.abs(w)<tol))
    nz=w[np.abs(w)>=tol]; gap=np.min(np.abs(nz)) if len(nz) else 0
    print(f" (2) eigenvalues of S∞ (dim {len(w)}): min={w[0]:.3e} max={w[-1]:.3e}; near-kernel dim={ker} (tol {tol:.1e}); gap above kernel={gap:.3e}")
    print(f"     smallest 6 |eig|: {np.sort(np.abs(w))[:6]}")
    # (4) locality: bandedness (energy in |i-j|<=2) vs total; Toeplitz-ness
    absS=np.abs(Sinf); tot=np.sum(absS)
    band=sum(np.sum(np.abs(np.diag(Sinf,k))) for k in range(-2,3))
    print(f" (4) locality: band(|i-j|<=2) energy / total = {band/tot:.3f}  (→1 local/banded; small=dense/nonlocal)")
print("\nREAD: S∞ converges (1) with a fixed near-kernel (the radical). If gap above kernel (2) is O(1) and STABLE")
print("across c, and S∞ is LOCAL/recognizable (4) => a finite polarization plausibly exists (construction has legs).")
print("If gap→0 with c OR S∞ dense/nonlocal => the quotient is also knife-edge/nonlocal = the same wall on S∞.")

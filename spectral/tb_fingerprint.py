#!/usr/bin/env python3
"""T-b FINGERPRINT TEST (§ bridge, 2026-08-24): does the arch-prime Weil quadratic form Q have the
SAME signature fingerprint as our verified arithmetic Hodge form?

Hodge form (Lean §1710/§1713/§1720):  hodgeForm = deg⊗deg − pairing  = (rank-1 POSITIVE) − (posdef),
  signature (1, dim−1): exactly ONE positive direction (= the degree functional), rest negative.

Bridge claim to test on the analytic side: the Weil form Q = arch − prime (Suzuki exact assembly)
  should have signature (dim−1, 1): exactly ONE NEGATIVE direction, and that direction = the s=1 POLE
  (Qpole, rank≤2). Then −Q has signature (1, dim−1) = the Hodge shape, and the single special mode of
  each (Weil's 1 neg = the pole ; Hodge's 1 pos = the degree) are the SAME object. That identity is
  exactly the T-b fingerprint the geometric object over Spec ℤ must realize.

We MEASURE (verifiable, no assumptions):
  (1) signature (n+, n-, n0) of Q;
  (2) whether the negative content is rank-1 and its eigenvector lives in the pole subspace;
  (3) whether the archimedean-only form (pole+core+tail) is PSD (arsenal §964-967 said yes);
  (4) whether removing the pole (project out Qpole's range) makes Q positive semidefinite.
"""
import numpy as np, sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
from suzuki_exact_Q import build_Q

def sig(w, tol=1e-8):
    return int((w > tol).sum()), int((w < -tol).sum()), int((np.abs(w) <= tol).sum())

def subspace_overlap(vec, basis_cols):
    """fraction of ||vec||^2 lying in span(basis_cols); basis_cols: (K,r)."""
    Qb, _ = np.linalg.qr(basis_cols)
    proj = Qb @ (Qb.T @ vec)
    return float(proj @ proj) / float(vec @ vec)

print("="*88)
print("T-b FINGERPRINT: signature of the arch-prime Weil form Q vs the Hodge (1,dim-1) shape")
print("="*88)
for a in [1.0, 1.5, 2.0, 2.5, 3.0]:
    K = max(70, int(30*a) + 40)
    Q, P = build_Q(a, K, return_pieces=True)
    Q = (Q + Q.T) / 2
    w, V = np.linalg.eigh(Q)
    npos, nneg, nz = sig(w * a)  # scale-free tol via *a is minor; use raw
    npos, nneg, nz = sig(w)
    # archimedean-only form
    Qarch = (P["pole"] + P["core"] + P["tail"]); Qarch = (Qarch + Qarch.T) / 2
    wa = np.linalg.eigvalsh(Qarch)
    apos, aneg, az = sig(wa)
    # pole subspace (range of Qpole, rank<=2)
    pole = (P["pole"] + P["pole"].T) / 2
    wp, Vp = np.linalg.eigh(pole)
    pole_cols = Vp[:, np.abs(wp) > 1e-9]          # columns spanning pole range
    # overlap of each negative eigenvector of Q with the pole subspace
    negidx = np.where(w < -1e-8)[0]
    overlaps = [subspace_overlap(V[:, i], pole_cols) for i in negidx] if len(negidx) else []
    # project Q onto complement of pole range; is it PSD?
    if pole_cols.shape[1] > 0:
        Qc, _ = np.linalg.qr(pole_cols)
        Pperp = np.eye(K) - Qc @ Qc.T
        Qproj = Pperp @ Q @ Pperp
        wproj = np.linalg.eigvalsh((Qproj + Qproj.T) / 2)
        # ignore the ~rank(pole) near-zero modes we projected out
        wproj_clean = np.sort(wproj)[pole_cols.shape[1]:]
        proj_min = wproj_clean.min()
    else:
        proj_min = float("nan")
    print(f"\na={a}  K={K}")
    print(f"  Q signature (n+,n-,n0)      = ({npos},{nneg},{nz})   λ_min/a={w[0]/a:+.4e}  λ_max/a={w[-1]/a:+.3e}")
    print(f"  arch-only (pole+core+tail)  = ({apos},{aneg},{az})   [arsenal: arch alone PSD ⟹ expect n-=0]")
    print(f"  pole rank                   = {pole_cols.shape[1]}")
    print(f"  neg-eigvec overlap w/ pole  = {['%.4f'%o for o in overlaps]}   [expect ~1.0 if neg mode = pole]")
    print(f"  Q on pole^perp: λ_min/a     = {proj_min/a:+.4e}   [expect ≥0 if removing pole ⟹ PSD]")

print("\n" + "="*88)
print("READ: fingerprint MATCH iff  n- of Q = 1 (mostly),  neg mode overlaps pole ~1, and")
print("Q restricted to pole^perp is PSD (λ_min≥0). That = the Weil form is 'rank-1 pole-negative")
print("+ posdef', the sign-flip of Hodge's 'rank-1 degree-positive + negdef'. Same special mode.")
print("="*88)

#!/usr/bin/env python3
"""§457 (is the coupling exploitable?): §456 localized the RH crux to prime's OFF-DIAGONAL coupling of archimedean
modes. Extract the coupling matrix C = (prime in archpole eigenbasis) with diagonal zeroed. Check: rank, norm vs
the diagonal margin, and whether it's structured (banded = nearest-mode coupling, or low-rank). Low-rank/banded =
a proof handle (finite-rank perturbation of the positive diagonal); full/generic = the coupling is irreducible."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0); w=3.0
for s in (4.0,3.5,3.0):
    cs=L.probes(BAND,s); f=L.build_forms(cs,w)
    AP=f['archpole']; PR=f['prime']; lam,U=np.linalg.eigh(AP)
    Pp=U.T@PR@U; n=len(cs)
    C=Pp-np.diag(np.diag(Pp))                  # off-diagonal coupling
    # full in archpole basis = diag(lam) - Pp; diagonal part positive, C is the threat
    Dfull=lam-np.diag(Pp)                        # diagonal of full (all positive per §456)
    lmin_diag=Dfull.min()
    lmin_full=float(np.linalg.eigvalsh(np.diag(lam)-Pp)[0])
    # coupling structure: bandedness (mass near diagonal) and effective rank
    absC=np.abs(C)
    band1=sum(absC[i,i+1]+absC[i+1,i] for i in range(n-1))/absC.sum()   # nearest-neighbor fraction
    svd=np.linalg.svd(C,compute_uv=False)
    effrank=int(np.sum(svd>0.1*svd[0]))
    print(f"s={s} n={n}: diag(full)min={lmin_diag:+.3e} λmin(full)={lmin_full:+.3e} (coupling reduces by {lmin_diag-lmin_full:.3e})")
    print(f"        ||C||={svd[0]:.3e}, effrank(C)={effrank}/{n}, nearest-mode fraction={band1:.2f}, singvals={np.array2string(svd[:5],precision=3)}")
print("""
reading: if C is LOW-RANK (effrank≪n) or BANDED (nearest-mode fraction high), the coupling is structured ⟹ a proof
could treat it as a finite-rank/tridiagonal perturbation of the positive diagonal (a real handle). If effrank~n
and coupling is spread (low nearest-mode fraction), the coupling is FULL/generic ⟹ irreducible, no cheap handle,
and the global control (§451/422) is genuinely needed. The gap diag(full)min − λmin(full) = how much the coupling
erodes = the quantity a proof must bound.
""")

#!/usr/bin/env python3
"""ATTACK A10: is the high-freq block Hi PROVABLY positive-definite via a diagonal scaling?
A symmetric matrix with positive diagonal is PD if it's an H-matrix (generalized diagonally
dominant): exists diagonal D>0 with Hi*D diag-dominant. Test: comparison matrix M(Hi) has
M_ii=|Hi_ii|, M_ij=-|Hi_ij|. Hi is an H-matrix iff M(Hi) is a nonsingular M-matrix iff
x=M(Hi)^{-1} 1 > 0 (all positive). If yes -> Hi>=0 is PROVABLE by manifest criterion ->
the positivity wall reduces to the low-freq Schur complement S_inf alone. Real reduction."""
import numpy as np, sys
sys.path.insert(0,"operators")
from critical_eigvec_id import build_tau_np

for c in [14]:
  for N in [96,160,200]:
    for M in [11,13,15]:
        tau,ns=build_tau_np(N,c)
        hi=np.where(np.abs(ns)>M)[0]
        Hi=tau[np.ix_(hi,hi)]
        diag=np.abs(np.diag(Hi))
        Comp=-np.abs(Hi); np.fill_diagonal(Comp,diag)   # comparison matrix M(Hi)
        try:
            x=np.linalg.solve(Comp,np.ones(len(Hi)))
            hmat = bool(np.all(x>0))
            xmin=float(x.min())
        except Exception as e:
            hmat=False; xmin=float('nan')
        # also: is Hi already plain diag-dominant on any row?
        dd_rows = int(np.sum(diag >= (np.sum(np.abs(Hi),axis=1)-diag)))
        print(f"c={c} N={N} M={M}: Hi is {len(Hi)}x{len(Hi)}  H-matrix(provably PD)={hmat}  "
              f"min(M^-1 1)={xmin:.2e}  plain-diagdom rows={dd_rows}/{len(Hi)}")

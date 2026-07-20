#!/usr/bin/env python3
"""
low_freq_schur.py -- ATTACK A9 (spawned from A2): localize the positivity wall to a FIXED
finite low-frequency block via Schur complement.

A2 showed the near-kernel is 23-dim degenerate -> no single critical eigenvector. But the
near-kernel = the LOW-FREQUENCY band (modes |n|<=~11). Hypothesis: tau PSD decomposes as
  - HIGH-freq block Hi (|n|>M): archimedean-dominated, MANIFESTLY PD (diagonally dominant?),
  - LOW-freq block Lo (|n|<=M): the hard part,
and tau>=0  <=>  Hi>=0  AND  Schur complement  S = Lo - Lo_Hi Hi^{-1} Hi_Lo  >= 0.
If (1) Hi is manifestly PD for all N, and (2) S converges to a FIXED finite PD matrix S_inf
as N->inf, then proving RH reduces to proving ONE small fixed matrix S_inf >= 0. That would
cut the infinite-dim wall down to a finite object -- a real decomposition, not a stat.

Tests: for M in a few values, across N: is Hi PD w/ a real gap? diagonally dominant? is S PD?
does S converge entrywise as N grows?
"""
import numpy as np
import json, sys
sys.path.insert(0,"/Users/scones/projects/riemann/operators")
from critical_eigvec_id import build_tau_np

def analyze(N, c, M):
    tau,ns = build_tau_np(N,c)
    lo = np.where(np.abs(ns)<=M)[0]
    hi = np.where(np.abs(ns)>M)[0]
    Lo = tau[np.ix_(lo,lo)]; Hi = tau[np.ix_(hi,hi)]
    LoHi = tau[np.ix_(lo,hi)]
    # high-freq block properties
    evHi = np.linalg.eigvalsh(Hi)
    # diagonal dominance margin of Hi: min_i (|Hi_ii| - sum_{j!=i}|Hi_ij|)
    dd = np.min(np.abs(np.diag(Hi)) - (np.sum(np.abs(Hi),axis=1)-np.abs(np.diag(Hi))))
    # Schur complement of Hi in tau (acting on Lo space)
    S = Lo - LoHi @ np.linalg.solve(Hi, LoHi.T)
    evS = np.linalg.eigvalsh(S)
    return {"D":2*N+1,"dimLo":len(lo),"dimHi":len(hi),
            "Hi_min_eig":float(evHi.min()),"Hi_max_eig":float(evHi.max()),
            "Hi_diag_dominance_margin":float(dd),
            "S_min_eig":float(evS.min()),"S_max_eig":float(evS.max()),
            "S":S}

if __name__=="__main__":
    c=int(sys.argv[1]) if len(sys.argv)>1 else 14
    Ns=[int(x) for x in sys.argv[2].split(",")] if len(sys.argv)>2 else [64,128,200]
    Ms=[int(x) for x in sys.argv[3].split(",")] if len(sys.argv)>3 else [11,13,15]
    out={"c":c,"results":{}}
    for M in Ms:
        print(f"\n{'='*70}\nBAND M={M} (low block = {2*M+1} modes |n|<={M})\n{'='*70}")
        Smats={}
        for N in Ns:
            r=analyze(N,c,M); Smats[N]=r.pop("S")
            print(f" N={N:4d} (D={r['D']}): Hi[{r['dimHi']}] min_eig={r['Hi_min_eig']:.3e} "
                  f"diagdom_margin={r['Hi_diag_dominance_margin']:.3e} | "
                  f"Schur S[{r['dimLo']}] min_eig={r['S_min_eig']:.3e}")
            out["results"].setdefault(str(M),{})[str(N)]={k:v for k,v in r.items()}
        # convergence of S across N (entrywise, on the common 2M+1 block)
        Ns_ok=[N for N in Ns if N in Smats]
        if len(Ns_ok)>=2:
            diffs=[float(np.max(np.abs(Smats[Ns_ok[i+1]]-Smats[Ns_ok[i]]))) for i in range(len(Ns_ok)-1)]
            print(f"  S entrywise max-change between consecutive N: {[f'{d:.3e}' for d in diffs]}"
                  f"   {'CONVERGING' if diffs[-1]<diffs[0] else 'not converging'}")
            out["results"][str(M)]["S_convergence_maxdiff"]=diffs
    p="/Users/scones/projects/riemann/results/low_freq_schur.json"
    json.dump(out,open(p,"w"),indent=2, default=lambda o: None)
    print(f"\nsaved -> {p}")
    print("\nREAD: if Hi min_eig>0 with a real gap AND diagdom_margin>0 (manifest PD) AND")
    print("Schur S min_eig>0 AND S converges -> wall reduced to fixed finite S_inf>=0. Proof lead.")

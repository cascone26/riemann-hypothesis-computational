#!/usr/bin/env python3
"""§404 (handoff ii): is the rank-collapse at RvM UNCONDITIONAL (counting, no arithmetic)? And exactly
how much of the corner is 'sign' (needs arithmetic) vs 'density' (free)?
rank(B_zero) = rank(2 sum_gamma F(gamma)F(gamma)^T) <= #zeros in support = RvM (unconditional).
Test: (a) numerical rank of the ANALYTIC side (primes-only) tracks RvM exactly (density=free);
(b) the SIGN residual = does the form stay PSD once rank saturates -> that is the ONLY arithmetic part.
Decompose lmin drop into rank-loss (counting) + sign-margin (arithmetic) across the plunge."""
import numpy as np
import l2_edge_sweep as L
w = 3.0
def rvm(T): 
    return T/(2*np.pi)*np.log(T/(2*np.pi))-T/(2*np.pi)+7/8 if T>2*np.pi else 0
print("rank(analytic) vs RvM count, + PSD-margin-above-rank (the arithmetic residual). band [14,44] w=3.")
print(f"{'s':>5} {'nprobe':>6} {'rvm_eff':>7} {'rank(B)':>7} {'rank=RvM?':>9} {'lmin/||B||':>11} {'2nd-smallest':>12}")
for s in np.arange(6.0, 2.4, -0.3):
    cs = L.probes(L.BAND, s)
    f = L.build_forms(cs, w)
    ev = np.linalg.eigvalsh(f['full']); nb = np.linalg.norm(f['full'],2)
    rk = int((ev > 1e-9*nb).sum())
    rvmeff = rvm(44+2*w) - rvm(14-2*w)
    print(f"{s:5.1f} {len(cs):6d} {rvmeff:7.2f} {rk:7d} {str(abs(rk-round(rvmeff))<=1):>9} "
          f"{ev[0]/nb:+11.2e} {ev[1]/nb:+12.2e}", flush=True)
print("\nreading: rank(B) tracks RvM (density=unconditional counting); the smallest eigenvalue staying")
print(">= 0 AS rank saturates is the ONLY place arithmetic (Weil positivity) enters = the sign half.")
print("If lmin>=0 exactly at rank saturation (not below), the corner arithmetic = 'PSD survives to the")
print("rank edge' — a codim-1 sign condition, everything else is free counting.")

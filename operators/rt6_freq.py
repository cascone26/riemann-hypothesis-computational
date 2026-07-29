#!/usr/bin/env python3
"""§417C: is the near-orthogonality (§417B) a PROVABLE oscillation mechanism? Hypothesis: v* (archimedean-
critical dir) is SMOOTH/low-frequency in r; the STRONG prime modes are HIGH-frequency (they carry cos(r log p)
for large p). Smooth × oscillatory ≈ 0 by Riemann-Lebesgue — UNCONDITIONAL cancellation. Measure the
'frequency' of v* and of each prime eigenmode as the mean |d/dr| (discrete), i.e. how oscillatory the mode is
as a function of probe center r. If strong prime modes are the oscillatory ones and v* is the smoothest, the
orthogonality has an unconditional oscillation-cancellation core; the RESIDUAL (how far 'some cancellation'
falls short of 'exactly enough for RH') is the real gap."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0); w=3.0
def roughness(vec,cs):
    o=np.argsort(cs); v=vec[o]; c=cs[o]
    d=np.diff(v)/np.diff(c); return float(np.sqrt(np.mean(d**2)))/ (np.sqrt(np.mean(v**2))+1e-12)
for s in (4.5,4.0,3.6):
    cs=L.probes(BAND,s)
    f=L.build_forms(cs,w); AP=f['archpole']; PR=f['prime']
    va,vav=np.linalg.eigh(AP); vstar=vav[:,0]
    pv,pu=np.linalg.eigh(PR); order=np.argsort(pv)[::-1]
    print(f"\ns={s} n={len(cs)}: roughness(v*)={roughness(vstar,cs):.3f}  (higher=more oscillatory)")
    print(f"  {'rank':>4} {'primeλ':>8} {'roughness':>10} {'|⟨v*,mode⟩|':>11}")
    for r,i in enumerate(order[:6]):
        print(f"  {r+1:>4} {pv[i]:8.3f} {roughness(pu[:,i],cs):10.3f} {abs(float(vstar@pu[:,i])):11.3f}")
print("""
reading: if v* is the SMOOTHEST (lowest roughness) and the top-λ prime modes are the ROUGHEST (oscillatory),
  the small overlaps are Riemann-Lebesgue = a provable unconditional cancellation core. The gap to RH is then
  purely in the LOW-frequency prime modes (few of them) whose overlap with smooth v* is NOT killed by
  oscillation — a finite-dim residual = the genuine crux (ties to §416 coherent/low-freq fluctuation & §389).
""")

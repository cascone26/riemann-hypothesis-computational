#!/usr/bin/env python3
"""§429 (sharp new question): does the EDGE deficit ρ_edge=||prime||op/archmargin (measured at probe count =
actual zero count in the band) DECREASE toward 1 as the band height γ→∞? If ρ_edge→1, the unconditional op-norm
bound becomes asymptotically tight at high γ (a real quantitative result). If ρ_edge stays bounded away from 1,
the cancellation deficit persists. Uses the instrument's real zeros to define the edge per band."""
import numpy as np, l2_fejer as LF
ZER=LF.ZEROS if hasattr(LF,'ZEROS') else None
def nz(band): return int(np.sum((ZER>=band[0])&(ZER<=band[1])))
def rho_at_edge(band):
    N=nz(band)
    if N<3: return None
    # choose s so #probes ≈ N (the edge)
    best=None
    for s in np.linspace((band[1]-band[0])/(N+3),(band[1]-band[0])/max(1,(N-1)),40):
        cs=LF.probes(band,s)
        if len(cs)==N:
            f=LF.build(cs); am=float(np.linalg.eigvalsh(f['archpole'])[0])
            ev=np.linalg.eigvalsh(f['prime']); opn=max(abs(ev[0]),abs(ev[-1]))
            lmf=float(np.linalg.eigvalsh(f['full'])[0])
            best=(N,am,opn,opn/am if am>1e-15 else np.inf,lmf); break
    return best
print(f"§429 ρ_edge vs band height. zeros available up to {ZER.max():.1f}")
print(f"{'band':>14} {'γ_mid':>6} {'Nz':>3} {'archmargin':>11} {'||prime||op':>11} {'ρ_edge':>7} {'λmin(full)':>11}")
for band in [(14,44),(20,50),(30,60),(40,70),(50,80),(60,90),(70,99)]:
    r=rho_at_edge(band)
    if r: 
        N,am,opn,rho,lmf=r
        print(f"{str(band):>14} {(band[0]+band[1])/2:6.0f} {N:3d} {am:+11.4e} {opn:11.4e} {rho:7.3f} {lmf:+11.4e}",flush=True)
print("""
reading: ρ_edge decreasing toward 1 as γ_mid grows => the archimedean margin (∼log γ) catches up to the bounded
  prime op-norm; unconditional positivity becomes asymptotically tight at high γ — quantitative progress. ρ_edge
  flat/growing => the deficit persists; RH cancellation is irreducibly needed at every height.
""")

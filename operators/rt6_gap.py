#!/usr/bin/env python3
"""§415 (ATTACK, new target): QUANTIFY what the second space must supply. full=(arch+pole)−prime, RH⟺full≥0.
Split the positivity into (a) the UNCONDITIONAL part (Gershgorin/diagonal-dominance lower bound on λmin —
provable with no RH input) and (b) the OFF-DIAGONAL CANCELLATION the true form exploits to stay ≥0 (the
genuine arithmetic/RH content). diag_margin(M)=min_i(M_ii − Σ_{j≠i}|M_ij|) is the Gershgorin bound: if ≥0,
positivity is free; if <0 while λmin≥0, the form RELIES on off-diagonal cancellation by exactly |diag_margin|.
GAP = λmin − diag_margin = how much the true arithmetic buys you over the trivial bound. Track across the
sweep to the rank edge; the gap at the edge = the precise deficit the second space (Hodge index) must close."""
import numpy as np, l2_edge_sweep as L
def diagmargin(M):
    d=np.diag(M); off=np.sum(np.abs(M),1)-np.abs(d); return float(np.min(d-off))
BAND=(14.0,44.0)
print(f"§415 unconditional(Gershgorin) vs true positivity. band{BAND} w=3. GAP=how much off-diag arith buys.")
print(f"{'s':>5} {'n':>3} {'λmin(full)':>11} {'Gershgorin':>11} {'off-diag reliance':>17} {'ratio λ/|diag|':>14}")
for s in np.arange(6.5,2.9,-0.25):
    cs=L.probes(BAND,s)
    if len(cs)<4: continue
    f=L.build_forms(cs,3.0); M=f['full']
    lm=float(np.linalg.eigvalsh(M)[0]); dm=diagmargin(M)
    reliance = -dm if dm<0 else 0.0    # amount of off-diag cancellation RH needs
    ratio = lm/abs(dm) if abs(dm)>1e-15 else np.nan
    print(f"{s:5.1f} {len(cs):3d} {lm:+11.3e} {dm:+11.3e} {reliance:17.3e} {ratio:14.3f}",flush=True)
print("""
reading: Gershgorin(diag_margin) < 0 everywhere the primes bite => positivity is NEVER free; RH-on-this-band
  ALWAYS relies on off-diagonal arithmetic cancellation of size |diag_margin|. The ratio λmin/|diag_margin|
  measures how EXACTLY the cancellation lands: ratio→0 at the edge => the arithmetic cancellation becomes
  razor-tight (buys just enough to stay ≥0, no slack) exactly where rank saturates. THAT razor-tightness —
  cancellation that lands exactly at 0, never overshooting negative — is what the second space must PROVE.
""")

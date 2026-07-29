#!/usr/bin/env python3
"""§428 (quantitative frontier): map the UNCONDITIONAL-positivity region precisely. full=archpole−prime is
provably PSD (no RH) wherever ||prime||_op ≤ λmin(archpole) [archimedean positivity = C-C theorem; prime bounded
by its op-norm]. Sweep band location (low γ vs high γ) × density on the clean Fejér instrument. Find the ratio
ρ = ||prime||_op / λmin(archpole): ρ≤1 => positivity UNCONDITIONALLY proven there; ρ>1 => needs cancellation (RH).
Question: is there ANY regime (very low density / wide spacing / specific band) where ρ<1 (a genuine unconditional
positive zone), and how does the frontier scale? This is the honest 'how far unconditionally' boundary."""
import numpy as np, l2_fejer as LF
def stats(cs):
    f=LF.build(cs); AP=f['archpole']; PR=f['prime']
    am=float(np.linalg.eigvalsh(AP)[0])
    ev=np.linalg.eigvalsh(PR); opn=max(abs(ev[0]),abs(ev[-1]))
    lmf=float(np.linalg.eigvalsh(f['full'])[0])
    return am,opn,lmf
print("§428 unconditional-positivity frontier (Fejér). ρ=||prime||op/archmargin; ρ≤1 => PROVEN PSD w/o RH.")
print(f"{'band':>14} {'s':>5} {'n':>3} {'archmargin':>11} {'||prime||op':>11} {'ρ':>7} {'λmin(full)':>11} {'uncond?':>8}")
for band in [(14,44),(14,20),(30,44),(60,90),(6,14),(2,10)]:
    for s in (12.0,9.0,7.0,5.5,4.5,3.5):
        cs=LF.probes(band,s)
        if len(cs)<3: continue
        am,opn,lmf=stats(cs)
        rho=opn/am if am>1e-15 else np.inf
        print(f"{str(band):>14} {s:5.1f} {len(cs):3d} {am:+11.4e} {opn:11.4e} {rho:7.2f} {lmf:+11.4e} "
              f"{'YES' if rho<=1 else 'no':>8}",flush=True)
print("""
reading: any 'YES' => a genuine unconditional-positivity zone (||prime||op below the archimedean margin there).
  ρ growing with density/band => the frontier; the ρ→1 crossover is the precise edge of what's provable without
  RH via the operator-norm bound. If ρ>1 EVERYWHERE the primes are resolved, the op-norm bound alone never proves
  positivity — matching that RH genuinely needs the cancellation (C-C combination), and quantifying the deficit ρ.
""")

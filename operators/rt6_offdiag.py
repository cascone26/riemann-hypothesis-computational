#!/usr/bin/env python3
"""§415b: is the razor-tight off-diagonal cancellation (§415) the PRIME PAIR-CORRELATION? Decompose the
off-diagonal deficit of full=pole+arch−prime into its arch vs prime parts. If the GROWING deficit (toward
the edge) is prime-dominated, then |diag_margin| = the prime pair-correlation mass, and RH-positivity ⟺
the prime pair-correlation cancellation lands exactly tight = Montgomery pair-correlation in the form's basis."""
import numpy as np, l2_edge_sweep as L
def offmass(M):  # total off-diagonal L1 mass
    return float(np.sum(np.abs(M))-np.sum(np.abs(np.diag(M))))
BAND=(14.0,44.0)
print(f"§415b off-diagonal mass: arch vs prime. band{BAND} w=3.")
print(f"{'s':>5} {'n':>3} {'offmass(arch)':>13} {'offmass(prime)':>14} {'prime/arch':>11} {'|diagmargin|':>12}")
for s in np.arange(6.0,2.9,-0.3):
    cs=L.probes(BAND,s)
    if len(cs)<4: continue
    f=L.build_forms(cs,3.0)
    pole=f['archpole']-((f['archpole'])-  # careful: recompute pole then arch
        (2.0*np.real(np.outer(
            np.exp(-((0.5j-cs)**2)/18)+np.exp(-((0.5j+cs)**2)/18),
            np.exp(-((0.5j-cs)**2)/18)+np.exp(-((0.5j+cs)**2)/18)))))
    fp=np.exp(-((0.5j-cs)**2)/18)+np.exp(-((0.5j+cs)**2)/18); pole=2.0*np.real(np.outer(fp,fp))
    arch=f['archpole']-pole; prime=f['prime']; M=f['full']
    d=np.diag(M); dm=float(np.min(d-(np.sum(np.abs(M),1)-np.abs(d))))
    oa,op=offmass(arch),offmass(prime)
    print(f"{s:5.1f} {len(cs):3d} {oa:13.3e} {op:14.3e} {op/oa:11.2f} {abs(dm):12.3e}",flush=True)
print("\nreading: prime/arch >> 1 and growing => the off-diagonal deficit is PRIME-dominated; the razor-tight")
print("cancellation of §415 is the prime pair-correlation. RH-positivity ⟺ prime pair-correlation lands tight.")

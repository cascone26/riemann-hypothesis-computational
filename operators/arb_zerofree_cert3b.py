#!/usr/bin/env python3
"""RIEMANN8 — zero-free certification of the SECOND §476g span (6 probes) at the open rung.
The 4-probe span {14,22,30,38} passed zero-free (arb_zerofree_cert3.out, banked §788).
This drives the same certify() on the 6-probe span {14,20,26,32,38,44} (§476g conditional
min pivot 1.63) — if it passes, BOTH §476g PD spans are unconditional by the direct route.
Same machinery, no new mathematics: arb_zerofree_cert3.certify(1.32, centers, RC=400).
(BAND1 = 44 in arb_fejer_cert2 covers the largest center; tail bound valid.)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arb_zerofree_cert3 import certify

if __name__ == "__main__":
    ok = certify(1.32, [14, 20, 26, 32, 38, 44])
    if ok:
        print("\n*** SECOND SPAN ZERO-FREE: both §476g PD spans now unconditional by the direct route. ***")
    else:
        print("\nNOT CERTIFIED on 6-probe span at RC=400 — next: raise RC / prec, or strip-lemma route (1b″).")

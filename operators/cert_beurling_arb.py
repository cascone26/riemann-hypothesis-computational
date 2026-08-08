#!/usr/bin/env python3
"""§830: arb-certify the §829 Panel-1 result → HANDLE. Below-capacity Weil positivity is ROBUST
across the t* crossing: certify Q=pole+arch−prime PD (pencil min > 0) on the FIXED K=5 span
[14,21.5,29,36.5,44] at each 4W in the Panel-1 sweep. Reuses the §823/§788 zero-free arb machinery
(module state set — the §806 landmine). If the pencil-min certifies positive with margin far above
the arb radii, the float 0.15–0.19 curve is upgraded from REPORT to HANDLE."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import cert_subtstar as cs   # certQ(fourW, CEN, rc) — sets ac.W/ac.TW/ac.c_tri correctly
ctx.prec = 250

SPAN = [14.0, 21.5, 29.0, 36.5, 44.0]   # the §829 Panel-1 fixed K=5 span (below capacity)
print("=" * 92)
print("§830 arb-CERTIFY below-capacity Weil positivity across the t* crossing (K=5 fixed span)")
print("span:", SPAN, "  (float §829 margins: 0.151 → 0.191 → 0.188 → 0.184)")
print("=" * 92)
print(f"{'4W':>6} {'t*':>9} {'cert pencil-min':>18} {'maxrad':>10}  verdict")
for fourW in [1.30, 1.60, 1.90, 2.05]:
    ts = 2 * np.pi * np.exp(fourW)
    m, mr = cs.certQ(fourW, SPAN, rc=400)
    if m is not None:
        verdict = "CERTIFIED PD (Handle)"
        val = f"{m:.6f}"
    else:
        verdict = "radii-swallowed (undecidable)"
        val = "—"
    print(f"{fourW:>6.2f} {ts:>9.2f} {val:>18} {mr:>10.1e}  {verdict}", flush=True)
print("=" * 92)
print("If all rows CERTIFIED: below-capacity Weil positivity across the full t* crossing is a HANDLE —")
print("the margin does NOT thin at t*; §829 Panel 1 upgraded REPORT→HANDLE.")

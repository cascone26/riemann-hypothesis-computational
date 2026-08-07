#!/usr/bin/env python3
"""THE MEASURED MARGIN CURVE m(W) — Euler's instrument (roundtable §788 → banks as §790).

Run the zero-free direct certificate (pole + arch - prime, interval LDL^T, NO zero knowledge)
at a ladder of support widths 4W, fixed span {14,22,30,38}, RC=400, and chart min-pivot vs 4W.

What each outcome buys (all three valuable):
  (i)  each PASS extends the unconditional certified territory rung by rung — including across
       the log4 = 1.386 entry (first PRIME-POWER rung, weight L(4)/sqrt(4) = ln2/2 ~ 0.347,
       barely half of w3 ~ 0.634; the ladder's rungs are unequal and gapped — no rung at log 6);
  (ii) the empirical m(W) curve vs the §485c margin-collapse rate e^{-9.78*4W} — measuring the
       functor W -> C_W before proving things about it;
  (iii) the 4W where RC=400 direct certification FIRST FAILS is a measured coordinate: the exact
       frontier where the strip lemma + frame architecture must take over from finite computation.

Guards (§788 Red audit) are enforced inside certify() itself: centers <= BAND1=44, RC > 2*BAND1,
4W < log(NM+1). A FAIL here is a RESULT (the frontier coordinate), not an error.
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arb_zerofree_cert3 import certify

SPAN = [14, 22, 30, 38]
LADDER = [1.32, 1.35, 1.40, 1.45, 1.50, 1.60, 1.70, 1.85, 2.00]

results = []
for fourW in LADDER:
    ok = certify(fourW, SPAN)
    results.append((fourW, ok))
    print(f"  == ladder point 4W={fourW}: {'CERTIFIED' if ok else 'NOT CERTIFIED (frontier coordinate)'} ==",
          flush=True)

print("\n===== MARGIN CURVE m(4W), span {14,22,30,38}, RC=400, zero-free =====")
print("(read min-pivot lower bounds from the per-run lines above; log-slope vs the")
print(" §485c collapse rate e^(-9.78*4W) is the comparison that matters)")
for fourW, ok in results:
    marks = []
    if fourW >= math.log(3): marks.append("p3 in")
    if fourW >= math.log(4): marks.append("4=2^2 in")
    if fourW >= math.log(5): marks.append("p5 in")
    if fourW >= math.log(7): marks.append("p7 in")
    print(f"  4W={fourW:.2f}  {'PASS' if ok else 'FAIL'}   [{', '.join(marks)}]")

#!/usr/bin/env python3
"""m(W) curve, DENSE variant — Nyquist-matched spans (companion to rung_ladder_scan.py; §790).

The coarse-span scan came back nearly FLAT (2.024 -> 2.016 across the log4 entry): a fixed
well-separated 4-probe span does not see the §485c collapse — that lives in the dense direction
(§476f corner / §769 rank-deficiency / §774 Nyquist edge). This scan holds the WINDOW fixed
([14,44]) and DENSIFIES the span with W to track the collapsing margin: spacing chosen at a fixed
multiple of the Nyquist mean zero spacing 2*pi/log(T/2pi) at the window center T=29 (~3.7), scaled
down as the certified §775 edge law predicts trouble. Two families per 4W:
  N=6  (spacing 6.0 — §788-proven safe at 1.32; watch it move)
  N=9  (spacing 3.75 — just above the §774 edge; expected to approach failure)
FAIL/near-singular here is the MEASUREMENT (the §485c collapse made zero-free-rigorous), not a bug.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arb_zerofree_cert3 import certify_lambda

SPANS = {
    "N6 (sp 6.0)":  [14, 20, 26, 32, 38, 44],
    "N9 (sp 3.75)": [14, 17.75, 21.5, 25.25, 29, 32.75, 36.5, 40.25, 44],
}
LADDER = [1.32, 1.45, 1.60, 1.80, 2.00]

for fourW in LADDER:
    for label, span in SPANS.items():
        print(f"\n===== 4W={fourW}  span {label} =====", flush=True)
        try:
            ok, lam = certify_lambda(fourW, span)
            print(f"  == dense point 4W={fourW} {label}: "
                  f"{f'CERTIFIED lambda_min_lb={lam:.6e}' if ok else 'NOT CERTIFIED (collapse coordinate)'} ==", flush=True)
        except AssertionError as e:
            print(f"  == dense point 4W={fourW} {label}: GUARD REFUSED ({e}) ==", flush=True)

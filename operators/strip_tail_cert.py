#!/usr/bin/env python3
"""Arb certification of the strip-envelope tail lemma applied to §476c (RUNG_STRIP_LEMMA.md → §791).

Recomputes the §476c enclosure with the tail bound valid for ARBITRARY zero configurations in the
critical strip above G=800 (two-sided, cosh^4(W*sigma_max) envelope from the Fejer closed form),
instead of the original one-sided on-line tail. If the certified enclosure of Q_no3(v) stays < 0,
the load-bearing indefiniteness holds with hypothesis set: {zeros <= 800 on-line (PT), classical
counting bound} — nothing about zeros above 800 at all.

Reuses arb_rung_certA2's computation by import (it runs top-level: witness reconstruction,
closed-form prime-3 term, zero side over 491 zeros, |v|_1); then swaps the tail.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx

import arb_rung_certA2 as ca   # runs the full §476c computation on import (prints its own log)

ctx.prec = 200
W = ca.Wa
G = ca.G
sv = ca.sv
Gm = arb(G - 44)

# on-line tail magnitude (certA2's formula, one-sided there)
t_online = sv*sv * 4/(W**4) * 2*(arb(G).log() + arb(1)/2) / (3*Gm**3)
# strip envelope factor: cosh^4(W * 1/2)
strip_factor = (W/2).cosh()**4
t_strip = t_online * strip_factor
term_tail_strip = 2 * arb(0).union(t_strip).union(-t_strip)   # two-sided, per-lemma

Q_strip = ca.term_zero + term_tail_strip + ca.term_prime3

print("\n===== STRIP-ENVELOPE TAIL CERT (§791) =====")
print(f"cosh^4(W/2) strip factor: {strip_factor}")
print(f"on-line tail bound (one-sided, §476c): +{t_online}")
print(f"strip tail bound (two-sided): +/- {t_strip}")
print(f"Q_no3(v) with ARBITRARY strip configurations above G={G}: {Q_strip}")
print(f"explicit bounds: Q_strip in [{float(Q_strip.lower()):.4e}, {float(Q_strip.upper()):.4e}]")
if Q_strip < 0:
    print("CERTIFIED: Q_no3(v) < 0 with NO on-line hypothesis above 800.")
    print("Hypothesis set, complete: zeros<=800 on-line (Platt-Trudgian), N(t+1)-N(t)<=2 ln t. Nothing else.")
else:
    print("NOT certified under arbitrary strip configurations at G=800 — fall back to citing PT to 3e12")
    print("(tail ~1e-36, trivial). Result banks either way; the hypothesis set is enumerated in both.")

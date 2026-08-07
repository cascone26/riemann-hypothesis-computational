#!/usr/bin/env python3
"""G-raise of §791 (Gauss's demand, §791-addendum): recompute the strip-unconditional Q_no3(v)
enclosure with G=2000 (1517 zeros, completeness N(2000)~1517 checked) instead of G=800.
Tail shrinks ~(756/1956)^3; margin predicted to return to ~-8e-5 with room behind it.
Witness + prime term reused from certA2 (import recomputes them rigorously); zero side rebuilt
here over the 33-digit zeros_100k.txt list (ball radius 1e-20, conservative for 33 printed digits).
Hypothesis set of the result: {zero list COMPLETE to 2000 and all on-line (Platt-Trudgian),
N(t+1)-N(t) <= 2 ln t}. Nothing else.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx

import arb_rung_certA2 as ca   # witness va, probes cs, |v|_1 sv, prime-3 term (runs its own G=800 log)

ctx.prec = 200
W = ca.Wa; cs = ca.cs; va = ca.va; n = ca.n; sv = ca.sv
G = 2000

zpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "zeros_100k.txt")
zeros = []
with open(zpath) as f:
    for line in f:
        s = line.strip()
        if not s: continue
        g = arb(s, 1e-20)
        if float(g.mid()) > G: break
        zeros.append(g)
print(f"\nG-raise: {len(zeros)} zeros <= {G} loaded (expect 1517)")
assert len(zeros) == 1517, f"completeness mismatch: {len(zeros)} vs N(2000)~1517"

zero_sum = arb(0)
for gb in zeros:
    psi = arb(0)
    for j in range(n):
        psi = psi + va[j]*ca.fej_arb(cs[j], gb)
    zero_sum = zero_sum + psi*psi
term_zero = 2*zero_sum
print(f"zero side (gamma<={G}, {len(zeros)} zeros): {term_zero}")

Gm = arb(G-44)
t_online = sv*sv * 4/(W**4) * 2*(arb(G).log() + arb(1)/2) / (3*Gm**3)
strip_factor = (W/2).cosh()**4
t_strip = t_online * strip_factor
term_tail_strip = 2 * arb(0).union(t_strip).union(-t_strip)

Q_strip = term_zero + term_tail_strip + ca.term_prime3
print(f"\n===== §791 G-RAISE (G={G}) =====")
print(f"strip tail bound (two-sided): +/- {float(2*t_strip):.4e}   (was +/- 8.04e-5 at G=800)")
print(f"Q_no3(v), arbitrary strip configs above {G}: explicit [{float(Q_strip.lower()):.4e}, {float(Q_strip.upper()):.4e}]")
if Q_strip < 0:
    print("CERTIFIED: Q_no3(v) < 0. Hypothesis set: {zero list complete to 2000 + on-line (PT), counting bound}.")
else:
    print("NOT certified at G=2000 — investigate before banking.")

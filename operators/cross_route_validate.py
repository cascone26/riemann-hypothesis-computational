#!/usr/bin/env python3
"""Cross-route validation of the Weil Gram (Red's demand, roundtable §788 → banks as §789).

The SAME 6x6 Weil-form Gram at 4W=1.32, centers {14,20,26,32,38,44}, computed via two logically
independent routes that share NO integrand code path for the transcendental part:
  route P (pole side, zero-free): pole_entry + arch_entry_v3 - prime_entry   [arb_zerofree_cert3]
  route Z (zero side, identity):  M[j][k] = 2 sum_{gamma>0} F_j(gamma) F_k(gamma)  [Guinand-Weil]
      truncated at G=800 (491 zeros, ball r=1e-9, on-line by Platt-Trudgian to 3e12) + rigorous
      entry-wise tail: |F_j(g)F_k(g)| <= 4/(W^4 (g-44)^4) for g > G, density N(t+1)-N(t) <= 2 ln t
      => tail_bound = 2 * 4/W^4 * 2(ln G + 1/2)/(3 (G-44)^3)   [same derivation as certA2].

WHY: pivot agreement (2.024 vs 2.02, 1.631 vs 1.63) across the two routes is the one validation
interval arithmetic cannot give — a wrong integrand/normalization/weight gets certified to 55 digits
just as happily. Entry-by-entry agreement of the two Grams within their rigorous enclosures validates
the explicit-formula BOOKKEEPING as implemented.

SCOPE (honest): route Z's tail one-sidedness/positivity uses on-line-ness of zeros above G (PT to
3e12; assumption above). This is a VALIDATION of implementation, not a certificate — route P remains
the unconditional certificate. Tail applied two-sided (±) here for safety anyway.

PASS criterion: for every (j,k), the route-P interval and route-Z interval INTERSECT.
Also reports max |mid_P - mid_Z| and the tightest slack.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

fourW = 1.32
CENTERS = [14, 20, 26, 32, 38, 44]
RC = 400
G = 800

ctx.prec = 220
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
Wa = ac.W
cs = [arb(c) for c in CENTERS]
n = len(cs)

# ---------- route P: pole + arch - prime (zero-free) ----------
print(f"route P (pole side): {n}x{n} Gram at 4W={fourW}, RC={RC} ...", flush=True)
MP = [[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j, n):
        e = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(RC)) - ac.prime_entry(cs[j], cs[k])
        MP[j][k] = MP[k][j] = e
    print(f"  P row {j} done (diag rad {float(MP[j][j].rad()):.2e})", flush=True)

# ---------- route Z: 2 sum_gamma F_j(gamma) F_k(gamma) + tail ----------
print(f"route Z (zero side): 491 zeros <= {G} + rigorous tail ...", flush=True)
ZC = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros500.npy'))

def fej_arb(c, r):
    out = arb(0)
    for cc in (c, -c):
        y = Wa*(r-cc)
        if y.contains(arb(0)):
            out = out + arb(1).union((y.sin()/y)**2 if not (y.mid() == 0 and y.rad() == 0) else arb(1))
        else:
            out = out + (y.sin()/y)**2
    return out

MZ = [[arb(0)]*n for _ in range(n)]
nz = 0
for g in ZC:
    if g > G: break
    gb = arb(float(g), 1e-9)
    F = [fej_arb(cs[j], gb) for j in range(n)]
    for j in range(n):
        for k in range(j, n):
            MZ[j][k] = MZ[j][k] + F[j]*F[k]
    nz += 1
tail = 2 * (arb(4)/Wa**4) * 2*(arb(G).log() + arb(1)/2) / (3*arb(G-44)**3)
tail_pm = arb(0).union(tail).union(-tail)
for j in range(n):
    for k in range(j, n):
        MZ[j][k] = 2*MZ[j][k] + tail_pm
        MZ[k][j] = MZ[j][k]
print(f"  Z done ({nz} zeros; entry tail bound +/- {float(tail):.3e})", flush=True)

# ---------- compare ----------
print("\nentry-by-entry comparison (PASS = intervals intersect):")
all_ok = True; maxd = 0.0; minslack = float('inf')
for j in range(n):
    for k in range(j, n):
        P, Z = MP[j][k], MZ[j][k]
        d = abs(float(P.mid()) - float(Z.mid()))
        budget = float(P.rad()) + float(Z.rad())
        ok = bool(P.overlaps(Z))
        all_ok &= ok
        maxd = max(maxd, d)
        minslack = min(minslack, budget - d)
        print(f"  [{CENTERS[j]},{CENTERS[k]}] P={float(P.mid()):+.6e} Z={float(Z.mid()):+.6e} "
              f"|dP-dZ|={d:.2e} budget={budget:.2e} {'OK' if ok else 'FAIL'}")
print(f"\nmax mid-discrepancy: {maxd:.3e}   tightest slack (budget - discrepancy): {minslack:.3e}")
if all_ok:
    print("CROSS-ROUTE VALIDATION PASSED: all 21 independent entries agree within rigorous enclosures.")
    print("The explicit-formula bookkeeping (pole/arch/prime vs zero side) is validated as implemented.")
else:
    print("MISMATCH FOUND — do NOT bank; investigate the disagreeing entries before any further certs.")

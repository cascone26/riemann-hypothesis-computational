#!/usr/bin/env python3
"""CERTIFY X — the corner↔complement cross-term block (Riemann round-3 item 1a → banks as §797).

Split basis at 4W=1.32, t* = 2pi e^{1.32} ~ 23.5:
  corner     A-block: centers {14, 20}      (below t*, the undersampled/deficiency region)
  complement C-block: centers {26, 32, 38, 44} (above t*, frame-controlled)
Full zero-free Weil Gram Q (pole + arch − prime, §788 route) on the joint 6-probe basis, interval
arithmetic; then certify the three block quantities and the composition inequality:
  m = certified λ_min(A) ·  c = certified λ_min(C) ·  ‖X‖ ≤ certified Frobenius bound
  PD of the whole  ⇐  m·c > ‖X‖²   (2×2 block Schur/eigen bound for [[A,X],[Xᵀ,C]])
This is the composition's PILOT: same numbers the direct §788-addendum cert already certifies PD,
now re-derived through the split — the template for classes too large for one direct certificate.
Hypothesis set: NONE (zero-free route throughout).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

fourW = 1.32; RC = 400
CORNER = [14, 20]; COMPL = [26, 32, 38, 44]
ALL = CORNER + COMPL
ctx.prec = 220
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
cs = [arb(c) for c in ALL]; n = len(cs)
nA = len(CORNER)

print(f"joint {n}x{n} zero-free Gram at 4W={fourW} (split at t*~23.5: corner {CORNER} | compl {COMPL})...", flush=True)
Q = [[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j, n):
        e = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(RC)) - ac.prime_entry(cs[j], cs[k])
        Q[j][k] = Q[k][j] = e
    print(f"  row {j} done", flush=True)

def lam_min_lb(M, m):
    def ok(sh):
        d = [None]*m; L = [[arb(0)]*m for _ in range(m)]
        for i in range(m):
            sd = M[i][i] - sh
            for k in range(i): sd -= L[i][k]*L[i][k]*d[k]
            if not (sd > 0): return False
            d[i] = sd
            for j2 in range(i+1, m):
                v = M[j2][i]
                for k in range(i): v -= L[j2][k]*L[i][k]*d[k]
                L[j2][i] = v/d[i]
        return True
    if not ok(arb(0)): return None
    lo, hi = 0.0, min(float(M[i][i].mid()) for i in range(m))
    for _ in range(40):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    return lo

A = [[Q[i][j] for j in range(nA)] for i in range(nA)]
C = [[Q[i][j] for j in range(nA, n)] for i in range(nA, n)]
m_lb = lam_min_lb(A, nA)
c_lb = lam_min_lb(C, n - nA)
xf2 = arb(0)
for i in range(nA):
    for j in range(nA, n):
        xf2 += Q[i][j]**2
xF = float(xf2.sqrt().upper())

print(f"\n===== CROSS-TERM BLOCK CERT (§797) =====")
print(f"corner  A {CORNER}: certified λ_min(A) ≥ {m_lb:.6e}" if m_lb else "corner A NOT PD")
print(f"compl   C {COMPL}: certified λ_min(C) ≥ {c_lb:.6e}" if c_lb else "compl C NOT PD")
print(f"cross   X: certified ‖X‖_2 ≤ ‖X‖_F ≤ {xF:.6e}")
if m_lb and c_lb:
    lhs = m_lb * c_lb; rhs = xF*xF
    print(f"composition inequality: m·c = {lhs:.6e}  vs  ‖X‖² = {rhs:.6e}  → "
          f"{'HOLDS — split-route PD certified' if lhs > rhs else 'FAILS — cross-terms dominate at this split'}")
    if lhs > rhs:
        import math
        lam_joint = 0.5*((m_lb+c_lb) - math.sqrt((m_lb-c_lb)**2 + 4*rhs))
        print(f"implied joint λ_min ≥ {lam_joint:.6e}  (direct §795 measured 0.661 — gap = the split's price)")

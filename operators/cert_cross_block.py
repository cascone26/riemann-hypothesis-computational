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
X = [[Q[i][j] for j in range(nA, n)] for i in range(nA)]      # 2 x 4
m_lb = lam_min_lb(A, nA)
c_lb = lam_min_lb(C, n - nA)
xf2 = arb(0)
for i in range(nA):
    for j in range(n - nA):
        xf2 += X[i][j]**2
xF = float(xf2.sqrt().upper())
# EXACT ||X||_2: X is 2x4 -> XX^T is 2x2 symmetric; lam_max closed form (interval-rigorous)
XXt = [[sum((X[i][k]*X[j][k] for k in range(n - nA)), arb(0)) for j in range(nA)] for i in range(nA)]
a_, b_, c_ = XXt[0][0], XXt[0][1], XXt[1][1]
lam_max_XXt = (a_+c_)/2 + (((a_-c_)/2)**2 + b_**2).sqrt()
x2 = float(lam_max_XXt.upper())               # = ||X||_2^2, exact up to interval radii

print(f"\n===== CROSS-TERM BLOCK CERT (§797) =====")
print(f"corner  A {CORNER}: certified λ_min(A) ≥ {m_lb:.6e}" if m_lb else "corner A NOT PD")
print(f"compl   C {COMPL}: certified λ_min(C) ≥ {c_lb:.6e}" if c_lb else "compl C NOT PD")
print(f"cross   X: ‖X‖_F ≤ {xF:.6e} (lossy) | EXACT ‖X‖₂² = λmax(XXᵀ) ≤ {x2:.6e}")
if m_lb and c_lb:
    lhs = m_lb * c_lb
    print(f"product inequality: m·c = {lhs:.6e}  vs  ‖X‖₂² = {x2:.6e}  → "
          f"{'HOLDS' if lhs > x2 else 'FAILS (product test lossy — see Schur)'}")
# LOSSLESS: Schur complement A - X C^{-1} X^T, interval route:
# solve C Y = X^T column-wise via interval LDL^T solve, then Schur = A - X Y; PD(C) ∧ PD(Schur) ⟺ PD(joint)
nC = n - nA
dC = [None]*nC; LC = [[arb(0)]*nC for _ in range(nC)]
okC = True
for i in range(nC):
    sd = C[i][i]
    for k in range(i): sd -= LC[i][k]*LC[i][k]*dC[k]
    if not (sd > 0): okC = False; break
    dC[i] = sd
    for j2 in range(i+1, nC):
        v = C[j2][i]
        for k in range(i): v -= LC[j2][k]*LC[i][k]*dC[k]
        LC[j2][i] = v/dC[i]
if okC:
    def solveC(rhs):
        y = [None]*nC
        for i in range(nC):
            s = rhs[i]
            for k in range(i): s -= LC[i][k]*y[k]
            y[i] = s
        for i in range(nC): y[i] = y[i]/dC[i]
        z = [None]*nC
        for i in reversed(range(nC)):
            s = y[i]
            for k in range(i+1, nC): s -= LC[k][i]*z[k]
            z[i] = s
        return z
    Y = [solveC([X[i][j] for j in range(nC)]) for i in range(nA)]   # rows: C^{-1} X^T columns
    Schur = [[A[i][j] - sum((X[i][k]*Y[j][k] for k in range(nC)), arb(0)) for j in range(nA)] for i in range(nA)]
    s_lb = lam_min_lb(Schur, nA)
    if s_lb:
        print(f"SCHUR (lossless): λ_min(A − XC⁻¹Xᵀ) ≥ {s_lb:.6e} > 0 ∧ C PD → JOINT PD CERTIFIED via the split")
        print(f"(the composition works losslessly; the m·c > ‖X‖² product form is what loses it)")
    else:
        print("SCHUR complement not certified PD — split genuinely fails beyond product-form loss")

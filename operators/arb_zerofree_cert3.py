#!/usr/bin/env python3
"""Zero-free rigorous certification v3 — tightened arch tail, conditioned spans, parameterized W.
GOAL: certify pole + arch - prime PD on a finite Fejer span WITHOUT any zero knowledge (unconditional
finite theorem, CC-flavor), at:
  (i)  4W=0.95  (p=2 rung, gate — CC-proven regime)
  (ii) 4W=1.32  (PAST log3: primes {2,3} — an unconditional certified instance INSIDE the open rung)
Tail fix vs arb_fejer_cert2: for r > RC > 2*BAND1, |F_j(r)F_k(r)| <= 4/(W(r-B1))^4 (fej(y)<=1/y^2 per
mirror term, both mirrors bounded by the near one with factor 2 each), and |Omega(r)| <= ln(r)+2 there
(digamma asymptotic |Re psi(1/4+ir/2)| <= ln(r) + 2 for r >= 10, safe). So
  |tail| <= (1/2pi) * int_RC^inf 4/(W(r-B1))^4 (ln r + 2) dr
         <= (2/pi) * (ln RC + 2 + 1/3) / (3 W^4 (RC-B1)^3)   [ln r <= ln RC + (r-RC)/RC, integrated safely:
            use ln r <= ln RC * r/RC crude for r>RC>e — simpler: ln r + 2 <= (ln RC + 2) * (r-B1)/(RC-B1)
            for r >= RC (since RHS grows linearly, LHS logarithmically, and equality at r=RC holds up to
            (RC-B1)/(RC-B1)=1 factor) => integral <= (ln RC + 2)/(RC-B1) * int (r-B1)^{-3} dr
            = (ln RC + 2) / (2 (RC-B1)^3) ] -- coded with this safe linear-domination form.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac

PREC = 220
ctx.prec = PREC

def arch_entry_v3(cj, ck, RC):
    def integ(r, _):
        return ac.fej_acb(acb(cj), r)*ac.fej_acb(acb(ck), r)*ac.omega_acb(r)
    ctx.prec = 300
    # manual subdivision: short segments bound the integrator's complex rectangles away from the
    # digamma pole at Im r = 1/2 (tall rectangles at small W caused inf radii)
    I = acb(0)
    seg = arb(2); a = -RC
    while a < RC:
        b = a + seg
        if b > RC: b = RC
        I = I + acb.integral(integ, acb(a), acb(b))
        a = b
    ctx.prec = PREC
    assert I.imag.contains(arb(0))
    val = I.real/(2*arb.pi())
    W = ac.W; B1 = ac.BAND1
    # tightened tail (derivation in header): (1/2pi)*4/W^4*(ln RC+2)/(2 (RC-B1)^3), x2 mirrors safety
    tail = (arb(2)/arb.pi()) * (RC.log()+2) / (W**4 * 2 * (RC-B1)**3) * 2
    return val + arb(0).union(tail.union(-tail))

def certify(fourW, centers, RC=400):
    # GUARDS (§788 Red audit): the tail derivation in the header is only valid under these —
    # violating any of them would print CERTIFIED on an invalid certificate.
    assert max(centers) <= float(ac.BAND1), \
        f"centers must be <= BAND1={float(ac.BAND1)} (mirror-domination tail needs r-B1 > 0 margin); got {max(centers)}"
    assert RC > 2 * float(ac.BAND1), f"tail bound needs RC > 2*BAND1; got RC={RC}, BAND1={float(ac.BAND1)}"
    assert fourW < float(arb(ac.NM + 1).log()), \
        f"prime sum truncates at NM={ac.NM}: need 4W < log(NM+1)={float(arb(ac.NM+1).log()):.3f}; got {fourW}"
    ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
    RCa = arb(RC)
    cs = [arb(c) for c in centers]; n = len(cs)
    prims = [nn for nn in range(2, 20) if ac.LAM[nn] is not None and float(arb(nn).log()) < fourW]
    print(f"\n4W={fourW} (primes in support: {prims}), centers {centers}, RC={RC}:", flush=True)
    M = [[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j, n):
            e = ac.pole_entry(cs[j], cs[k]) + arch_entry_v3(cs[j], cs[k], RCa) - ac.prime_entry(cs[j], cs[k])
            M[j][k] = M[k][j] = e
        print(f"  row {j} done (diag rad {float(M[j][j].rad()):.2e})", flush=True)
    d = [None]*n; L = [[arb(0)]*n for _ in range(n)]
    for i in range(n):
        sd = M[i][i]
        for k in range(i): sd -= L[i][k]*L[i][k]*d[k]
        if not (sd > 0):
            print(f"  PIVOT {i}: {sd} not>0 -> NOT CERTIFIED"); return False
        d[i] = sd
        for j2 in range(i+1, n):
            v = M[j2][i]
            for k in range(i): v -= L[j2][k]*L[i][k]*d[k]
            L[j2][i] = v/d[i]
    mn = min(float(x.mid())-float(x.rad()) for x in d)
    print(f"  ALL {n} PIVOTS POSITIVE (min pivot lower bound {mn:.3e}) -> ZERO-FREE CERTIFIED PD", flush=True)
    return True

def _ldlt_all_pos(M, n, shift):
    d = [None]*n; L = [[arb(0)]*n for _ in range(n)]
    for i in range(n):
        sd = M[i][i] - shift
        for k in range(i): sd -= L[i][k]*L[i][k]*d[k]
        if not (sd > 0): return False
        d[i] = sd
        for j2 in range(i+1, n):
            v = M[j2][i]
            for k in range(i): v -= L[j2][k]*L[i][k]*d[k]
            L[j2][i] = v/d[i]
    return True

def certify_lambda(fourW, centers, RC=400):
    """Like certify() but returns (ok, certified λ_min lower bound) — Red's §788/§791 demand:
    pivots certify sign, not margin; the m(W) that enters any composition inequality must be λ_min.
    Bisection on shift: largest λ with interval-LDLᵀ(M−λI) all pivots > 0 is a CERTIFIED λ_min lb
    (if M−λI is PD then λ_min(M) > λ, rigorous in interval arithmetic)."""
    assert max(centers) <= float(ac.BAND1), f"centers must be <= BAND1={float(ac.BAND1)}"
    assert RC > 2 * float(ac.BAND1), f"need RC > 2*BAND1"
    assert fourW < float(arb(ac.NM + 1).log()), f"prime sum truncates at NM={ac.NM}"
    ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
    RCa = arb(RC)
    cs = [arb(c) for c in centers]; n = len(cs)
    prims = [nn for nn in range(2, 20) if ac.LAM[nn] is not None and float(arb(nn).log()) < fourW]
    print(f"\n4W={fourW} (primes in support: {prims}), centers {centers}, RC={RC}:", flush=True)
    M = [[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j, n):
            e = ac.pole_entry(cs[j], cs[k]) + arch_entry_v3(cs[j], cs[k], RCa) - ac.prime_entry(cs[j], cs[k])
            M[j][k] = M[k][j] = e
        print(f"  row {j} done (diag rad {float(M[j][j].rad()):.2e})", flush=True)
    if not _ldlt_all_pos(M, n, arb(0)):
        print("  NOT CERTIFIED PD (shift 0 fails)"); return False, None
    lo, hi = 0.0, min(float(M[i][i].mid()) for i in range(n))
    for _ in range(40):
        mid = (lo + hi)/2
        if _ldlt_all_pos(M, n, arb(mid)): lo = mid
        else: hi = mid
    print(f"  ZERO-FREE CERTIFIED PD, certified lambda_min lower bound {lo:.6e}", flush=True)
    return True, lo

if __name__ == "__main__":
    ok1 = certify(0.95, [14, 22, 30, 38])
    ok2 = certify(1.32, [14, 22, 30, 38])
    if ok2:
        print("\n*** UNCONDITIONAL CERTIFIED POSITIVITY INSTANCE INSIDE THE OPEN RUNG (4W>log3, primes {2,3},")
        print("    no zero knowledge anywhere) — a finite zero-free theorem on this span. ***")

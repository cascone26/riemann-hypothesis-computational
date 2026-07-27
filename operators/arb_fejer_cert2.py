#!/usr/bin/env python3
"""§382: RIGOROUS Fejer certificate, feasible via §381 closed-form prime side.
B_full = pole(closed) + arch(one acb.integral/entry, sinc^4*Omega + 1/r^4 tail) - prime(CLOSED, §381).
Only 36 arch integrals for n=8 (was 36x30 oscillatory in §380 -> infeasible). Interval LDL^T.
"""
import sys
from flint import arb, acb, ctx

PREC = 200
ctx.prec = PREC
W = arb(6) / 5           # 1.2
BAND0, BAND1 = arb(14), arb(44)
RC = arb(140)

def probes(s):
    cs = []; c = BAND0
    while c <= BAND1 + arb('1e-9'):
        cs.append(c); c = c + s
    return cs

def fej_acb(c, r):
    out = acb(0)
    for cc in (c, -c):
        y = W * (r - cc)
        sc = y.sinc()
        out = out + sc * sc
    return out

def omega_acb(r):
    ii = acb(0, 1); q = acb(arb(1) / 4)
    return ((q + ii * r / 2).digamma() + (q - ii * r / 2).digamma()) / 2 - arb.pi().log()

def arch_entry(cj, ck):
    def integ(r, _):
        return fej_acb(acb(cj), r) * fej_acb(acb(ck), r) * omega_acb(r)
    ctx.prec = 280
    I = acb.integral(integ, acb(-RC), acb(RC))
    ctx.prec = PREC
    assert I.imag.contains(arb(0)), "arch imag"
    val = I.real / (2 * arb.pi())
    A = 1 / (W * (RC - BAND1)) ** 2
    tail = 4 * (2 * A) * (2 * A) * (RC.log() + 2) * RC / (2 * arb.pi())
    return val + arb(0).union(tail.union(-tail))

def pole_entry(cj, ck):
    z = acb(0, arb(1) / 2)
    return 2 * (fej_acb(acb(cj), z) * fej_acb(acb(ck), z)).real

# ---- closed-form prime (arb port of §381) ----
TW = 2 * W; c_tri = arb.pi() / W
def Jexact_arb(u, d):
    lo = max(-TW, -TW - u); hi = min(TW, TW - u)
    if not (hi > lo): return acb(0)
    cand = [lo, hi]
    for x in (arb(0), -u):
        if lo < x < hi: cand.append(x)
    bpts = sorted(cand, key=lambda t: float(t.mid()))
    tot = acb(0)
    for a_, b_ in zip(bpts[:-1], bpts[1:]):
        mid = (a_ + b_) / 2
        sv = arb(1) if mid > 0 else arb(-1)
        A1 = c_tri; B1 = -c_tri * sv / TW
        su = arb(1) if (u + mid) > 0 else arb(-1)
        A2 = c_tri * (1 - su * u / TW); B2 = -c_tri * su / TW
        P0 = A1 * A2; P1 = A1 * B2 + B1 * A2; P2 = B1 * B2
        def aq(v):
            if d == 0 or (hasattr(d, 'is_zero') and False):
                pass
            e = acb(0, -1) * d      # -i d
            # handle d==0 via arb check
            return None
        # do d==0 vs not explicitly (d is arb)
        if d.contains(arb(0)) and d.mid() == 0 and d.rad() == 0:
            def anti(v): return acb(P0 * v + P1 * v * v / 2 + P2 * v * v * v / 3)
        else:
            def anti(v):
                ee = (acb(0, -1) * d * v).exp()
                r = acb(0, 1) * P2 / d
                q = (2 * r - P1) / (acb(0, 1) * d)
                p = (q - P0) / (acb(0, 1) * d)
                return ee * (p + q * v + r * v * v)
        tot = tot + anti(b_) - anti(a_)
    return tot

def vm():
    NM = 200; lam = [None] * (NM + 1); comp = [False] * (NM + 1)
    for p in range(2, NM + 1):
        if not comp[p]:
            for m in range(2 * p, NM + 1, p): comp[m] = True
            lp = arb(p).log(); pk = p
            while pk <= NM: lam[pk] = lp; pk *= p
    return lam, NM
LAM, NM = vm()

def gpair_closed_arb(cj, ck, u):
    tot = acb(0)
    for a in (cj, -cj):
        for b in (ck, -ck):
            tot = tot + (acb(0, 1) * u * b).exp() * Jexact_arb(u, a - b)
    return tot.real / (2 * arb.pi()) ** 2

def prime_entry(cj, ck):
    s = arb(0); fourw = 4 * W
    for n in range(2, NM + 1):
        if LAM[n] is None: continue
        u = arb(n).log()
        if u > fourw: continue
        s = s + LAM[n] / arb(n).sqrt() * gpair_closed_arb(cj, ck, u)
    return 2 * s

def certify(s_str):
    s = arb(s_str); cs = probes(s); n = len(cs)
    print(f"\ns={s_str}: {n} Fejer probes — rigorous B_full (closed prime + {n*(n+1)//2} arch integrals)...", flush=True)
    M = [[None] * n for _ in range(n)]
    for j in range(n):
        for k in range(j, n):
            M[j][k] = M[k][j] = pole_entry(cs[j], cs[k]) + arch_entry(cs[j], cs[k]) - prime_entry(cs[j], cs[k])
        print(f"  row {j} done", flush=True)
    d = [None] * n; Lm = [[arb(0)] * n for _ in range(n)]
    for i in range(n):
        sd = M[i][i]
        for k in range(i): sd -= Lm[i][k] * Lm[i][k] * d[k]
        if not (sd > 0):
            print(f"  PIVOT {i}: {sd} not >0 -> NOT CERTIFIED s={s_str}"); return False
        d[i] = sd
        for j2 in range(i + 1, n):
            v = M[j2][i]
            for k in range(i): v -= Lm[j2][k] * Lm[i][k] * d[k]
            Lm[j2][i] = v / d[i]
    print(f"  ALL {n} PIVOTS POSITIVE -> Fejer B_full(s={s_str}) CERTIFIED PD (2nd-instrument THEOREM)")
    print(f"  smallest pivot {min(d, key=lambda x: x.mid())}")
    return True

if __name__ == "__main__":
    for s_str in (sys.argv[1:] or ["4.0"]):
        certify(s_str)

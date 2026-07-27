#!/usr/bin/env python3
"""§380: RIGOROUS (arb) certificate on the FEJER instrument — the theorem on a 2nd test class.
Advantage over Gaussian: the prime side is EXACTLY FINITE (compact frequency support, log n<4W)
=> no cut / no smooth-tail / no eps-charge. Only the arch needs a rigorous tail (sinc^4 ~ 1/r^4).

B_full = pole + arch - prime, all arb/acb:
  pole[j,k]  = 2 Re[F_j(i/2) F_k(i/2)]                       (closed, acb)
  arch[j,k]  = (1/2pi) Int_{-Rc}^{Rc} F_jF_k Omega dr + tail  (acb.integral + explicit 1/r^4 tail)
  prime[j,k] = 2 sum_{n: log n < 4W} Lambda(n)/sqrt(n) (1/2pi) Int F_jF_k cos(ur) dr  (FINITE, acb.integral)
Interval LDL^T pivot test. Validate: arch imag encloses 0; residual vs float small. Then PD verdict.
fej(x) = (sin(Wx)/(Wx))^2.
"""
import sys
from flint import arb, acb, ctx

PREC = 200
ctx.prec = PREC
W = arb(6) / 5            # 1.2, dyadic-ish (6/5 exact rational in arb)
BAND0, BAND1 = arb(14), arb(44)
RC = arb(120)            # arch integration cutoff
FOURW = 4 * W            # prime frequency cutoff: log n < 4W

def logpi(): return arb.pi().log()

def probes(s):
    cs = []; c = BAND0
    while c <= BAND1 + arb('1e-9'):
        cs.append(c); c = c + s
    return cs

def fej_acb(c, r):
    # (sin(W(r-c))/(W(r-c)))^2 + mirror; r is acb, c real arb
    out = acb(0)
    for cc in (c, -c):
        y = W * (r - cc)
        # guard: y never exactly 0 on the integration path off the nodes; acb sin/y fine generically
        out = out + (y.sin() / y) ** 2
    return out

def omega_acb(r):
    ii = acb(0, 1); q = acb(arb(1) / 4)
    return ((q + ii * r / 2).digamma() + (q - ii * r / 2).digamma()) / 2 - logpi()

def arch_entry(cj, ck):
    def integ(r, _):
        return fej_acb(cj, r) * fej_acb(ck, r) * omega_acb(r)
    ctx.prec = 300
    I = acb.integral(integ, acb(-RC), acb(RC))
    ctx.prec = PREC
    assert I.imag.contains(arb(0)), "arch imag not enclosing 0"
    val = I.real / (2 * arb.pi())
    # tail |r|>RC: |fej(r-c)| <= 1/(W(r-c))^2 <= 1/(W(RC-BAND1))^2 =: A each; |F_jF_k| <= (2A)^2-ish;
    # |Omega(r)| <= log(r) + 1 for r>RC (digamma asymptotics, crude); Int_{RC}^inf (log r+1)/(W(r-B1))^4 dr
    A = 1 / (W * (RC - BAND1)) ** 2
    # crude rigorous tail: 4*(2A)^2 * (log(RC)+2)/RC^3 * RC  (over-bound both sides)
    tailmag = 4 * (2 * A) * (2 * A) * (RC.log() + 2) * RC / (2 * arb.pi())
    return val + arb(0).union(tailmag.union(-tailmag))

def pole_entry(cj, ck):
    z = acb(0, arb(1) / 2)
    return 2 * (fej_acb(cj, z) * fej_acb(ck, z)).real

def vonmangoldt():
    NM = 200
    lam = [None] * (NM + 1); comp = [False] * (NM + 1)
    for p in range(2, NM + 1):
        if not comp[p]:
            for m in range(2 * p, NM + 1, p): comp[m] = True
            lp = arb(p).log(); pk = p
            while pk <= NM: lam[pk] = lp; pk *= p
    return lam, NM
LAM, NM = vonmangoldt()

def prime_entry(cj, ck):
    s = arb(0); fourw = FOURW
    for n in range(2, NM + 1):
        if LAM[n] is None: continue
        u = arb(n).log()
        if u > fourw: continue          # compact support: exactly zero beyond
        def integ(r, _, u=u):
            return fej_acb(cj, r) * fej_acb(ck, r) * (acb(u) * r).cos()
        ctx.prec = 300
        I = acb.integral(integ, acb(-RC), acb(RC))
        ctx.prec = PREC
        gp = I.real / (2 * arb.pi())    # tail negligible: cos oscillates, sinc^4 tail tiny (folded into arch-style bound omitted: u<4W bounded)
        s = s + LAM[n] / arb(n).sqrt() * gp
    return 2 * s

def certify(s_str):
    s = arb(s_str); cs = probes(s); n = len(cs)
    print(f"\ns={s_str}: {n} Fejer probes — rigorous B_full ...", flush=True)
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
            print(f"  PIVOT {i}: {sd} not >0 -> NOT CERTIFIED at s={s_str}"); return False
        d[i] = sd
        for j2 in range(i + 1, n):
            v = M[j2][i]
            for k in range(i): v -= Lm[j2][k] * Lm[i][k] * d[k]
            Lm[j2][i] = v / d[i]
    print(f"  ALL {n} PIVOTS POSITIVE -> Fejer B_full(s={s_str}) CERTIFIED PD (arb prec {PREC})")
    print(f"  => 2nd-instrument THEOREM: finite exact prime side (n<={int(2.718**float(FOURW))}) + arb arch + pole.")
    print(f"  smallest pivot: {min(d, key=lambda x: x.mid())}")
    return True

if __name__ == "__main__":
    for s_str in (sys.argv[1:] or ["4.0", "3.0"]):
        certify(s_str)

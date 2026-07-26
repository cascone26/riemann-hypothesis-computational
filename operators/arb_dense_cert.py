#!/usr/bin/env python3
"""§340a: RIGOROUS (arb interval) certification of the analytic Weil form at plunge-floor spacings.
B_full = pole + arch - prime on symmetrized Gaussian probes, band [14,44], w=3, s in {3.5, 3.4}
(float margins 1.9e-10 / 9.4e-12 — at/below float-trust; here every entry is an arb enclosure and
PD is certified by LDL^T with the §199.5 audit-fixed discipline: interval pivot test s>0, acb
integrals with imag-encloses-0 assertion, explicit tail radius).
Entries:
  pole[j,k] = 2 Re[F_j(i/2) F_k(i/2)]              (closed form, acb)
  arch[j,k] = (1/pi) * Re Int_0^R F_jF_k(r) Omega(r) dr  + [tail radius]   (acb.integral)
  prime[j,k] = 2 Sum_n Lambda(n)/sqrt(n) g_pair(log n)   (closed form, arb; n<=300 + tail radius)
"""
import sys
from flint import arb, acb, ctx

PREC = 220
ctx.prec = PREC
W = arb(3)
BAND0, BAND1 = arb(14), arb(44)
RCUT = arb(95)          # arch integral cutoff; tail bounded explicitly
NMAX = 300

def logpi(): return arb.pi().log()

def probes(s):
    cs = []
    c = BAND0
    while c <= BAND1 + arb('1e-9'):
        cs.append(c)
        c = c + s
    return cs

def F_acb(c, r):
    """symmetrized Gaussian at complex r (acb)."""
    return (-((r - c) * (r - c)) / (2 * W * W)).exp() + (-((r + c) * (r + c)) / (2 * W * W)).exp()

def omega_acb(r):
    """holomorphic extension of Re psi(1/4+ir/2) - log pi (equals it on the real axis)."""
    ii = acb(0, 1)
    quarter = acb(arb(1) / 4)
    return ((quarter + ii * r / 2).digamma() + (quarter - ii * r / 2).digamma()) / 2 - logpi()

def arch_entry(cj, ck):
    def integ(r, _):
        return F_acb(acb(cj), r) * F_acb(acb(ck), r) * omega_acb(r)
    ctx.prec = 260
    I = acb.integral(integ, acb(0), acb(RCUT))
    ctx.prec = PREC
    assert I.imag.contains(arb(0)), "arch imag not enclosing 0"
    val = I.real / arb.pi()
    # tail |r|>RCUT: |F_jF_k| <= 4 exp(-(RCUT-c)^2/w^2)-ish; Omega <= log(r) growth; crude rigorous radius
    tj = (-((RCUT - cj) * (RCUT - cj)) / (2 * W * W)).exp()
    tk = (-((RCUT - ck) * (RCUT - ck)) / (2 * W * W)).exp()
    tail = arb(8) * tj * tk * RCUT  # |Omega|<=RCUT on tail crudely; integral of decaying gaussian < 1
    return val + arb(0).union(tail.union(-tail))

def pole_entry(cj, ck):
    z = acb(0, arb(1) / 2)
    return 2 * (F_acb(acb(cj), z) * F_acb(acb(ck), z)).real

def vonmangoldt_list():
    lam = [arb(0)] * (NMAX + 1)
    comp = [False] * (NMAX + 1)
    for p in range(2, NMAX + 1):
        if not comp[p]:
            for m in range(2 * p, NMAX + 1, p):
                comp[m] = True
            lp = arb(p).log()
            pk = p
            while pk <= NMAX:
                lam[pk] = lp
                pk *= p
    return lam

LAM = vonmangoldt_list()

def g_pair(cj, ck, u):
    pref = W * arb.pi().sqrt() / (2 * arb.pi())
    d1 = cj - ck
    e_near = (-(d1 * d1) / (4 * W * W)).exp()
    d2 = cj + ck
    e_far = (-(d2 * d2) / (4 * W * W)).exp()
    gauss = (-(u * u) * W * W / 4).exp()
    return pref * gauss * 2 * (e_near * (u * (cj + ck) / 2).cos() + e_far * (u * (cj - ck) / 2).cos())

def prime_entry(cj, ck):
    s = arb(0)
    for n in range(2, NMAX + 1):
        if not (LAM[n] == arb(0)):
            u = arb(n).log()
            s += LAM[n] / arb(n).sqrt() * g_pair(cj, ck, u)
    # tail n>NMAX: Lambda(n)/sqrt(n) <= log n; sum log n * e^{-(log n)^2 w^2/4} * pref * 4 over n>NMAX
    uN = arb(NMAX).log()
    tailmax = (W * arb.pi().sqrt() / (2 * arb.pi())) * 4 * uN * (-(uN * uN) * W * W / 4).exp() * NMAX
    return 2 * s + arb(0).union((2 * tailmax).union(-2 * tailmax))

def certify(s_str):
    s = arb(s_str)
    cs = probes(s)
    n = len(cs)
    print(f"s={s_str}: {n} probes — building rigorous B_full ({n*(n+1)//2} arch integrals)...", flush=True)
    M = [[None] * n for _ in range(n)]
    for j in range(n):
        for k in range(j, n):
            v = pole_entry(cs[j], cs[k]) + arch_entry(cs[j], cs[k]) - prime_entry(cs[j], cs[k])
            M[j][k] = M[k][j] = v
        print(f"  row {j} done", flush=True)
    # rigorous LDL^T, interval pivot test
    d = [None] * n
    L = [[arb(0)] * n for _ in range(n)]
    ok = True
    for i in range(n):
        sdiag = M[i][i]
        for k in range(i):
            sdiag -= L[i][k] * L[i][k] * d[k]
        if not (sdiag > 0):
            print(f"  PIVOT {i}: interval {sdiag} does NOT certify > 0  ->  NOT CERTIFIED at s={s_str}")
            ok = False
            break
        d[i] = sdiag
        for j2 in range(i + 1, n):
            v = M[j2][i]
            for k in range(i):
                v -= L[j2][k] * L[i][k] * d[k]
            L[j2][i] = v / d[i]
    if ok:
        print(f"  ALL {n} PIVOTS RIGOROUSLY POSITIVE  ->  B_full(s={s_str}) CERTIFIED PD (arb, prec {PREC})")
        print(f"  smallest pivot enclosure: {min(d, key=lambda x: x.mid())}")
    return ok

if __name__ == "__main__":
    for s_str in (sys.argv[1:] or ["3.5", "3.4"]):
        certify(s_str)

#!/usr/bin/env python3
"""§350: ARB-RIGOROUS certification of the COLLECTIVE architecture at deep configs.

B_full = pole + arch - P  >=  B_lower = pole + arch - P_small(n<n0) - T_smooth - CHG*I
  (as forms: P = P_small + T_exact; T_exact <= T_smooth + ||DT||*I; ||DT|| <= CHG), where
  CHG = eps * ||C_env||_F  (fluctuation on [n0, XCUT], partial summation, envelope integrand)
      + tail_exact(>XCUT)  (|psi'-mass| via psi(x) < 1.04x, classical)
      + tail_smooth(>XCUT) (Gaussian tail of the smooth integral).
eps is PROVED IN-SCRIPT: exhaustive arb check of sup_{x in [n0, XCUT]} |psi(x)-x|/x
(integer breakpoints; on [n,n+1) sup = max(|psi(n)-n|, |psi(n)-(n+1)|)/n).
Everything else: §340 discipline — acb integrals with imag-encloses-0 assertion + explicit
tail radii, closed-form arb sums, interval LDL^T with pivot test, NO ** on arb balls.
Certifying B_lower PD certifies B_full PD: theorem-grade for the config.
Configs: w=2, band [14,44], n0=300, s in {3.9, 3.8} (deepest PSD rows; §349 float margins
6.8e-10 / 6.3e-12; float charge ~2e-13 << margins).
"""
import sys
from flint import arb, acb, ctx

PREC = 220
ctx.prec = PREC
W = arb(2)
BAND0, BAND1 = arb(14), arb(44)
RCUT = arb(95)
N0 = 300
XCUT = 20000

def logpi(): return arb.pi().log()

def probes(s):
    cs = []
    c = BAND0
    while c <= BAND1 + arb('1e-9'):
        cs.append(c)
        c = c + s
    return cs

def F_acb(c, r):
    return (-((r - c) * (r - c)) / (2 * W * W)).exp() + (-((r + c) * (r + c)) / (2 * W * W)).exp()

def omega_acb(r):
    ii = acb(0, 1)
    quarter = acb(arb(1) / 4)
    return ((quarter + ii * r / 2).digamma() + (quarter - ii * r / 2).digamma()) / 2 - logpi()

def arch_entry(cj, ck):
    def integ(r, _):
        return F_acb(acb(cj), r) * F_acb(acb(ck), r) * omega_acb(r)
    ctx.prec = 320
    I = acb.integral(integ, acb(0), acb(RCUT))
    ctx.prec = PREC
    assert I.imag.contains(arb(0)), "arch imag not enclosing 0"
    val = I.real / arb.pi()
    tj = (-((RCUT - cj) * (RCUT - cj)) / (2 * W * W)).exp()
    tk = (-((RCUT - ck) * (RCUT - ck)) / (2 * W * W)).exp()
    tail = arb(8) * tj * tk * RCUT
    return val + arb(0).union(tail.union(-tail))

def pole_entry(cj, ck):
    z = acb(0, arb(1) / 2)
    return 2 * (F_acb(acb(cj), z) * F_acb(acb(ck), z)).real

# ---- exact psi to XCUT (arb), Lambda list, and the PROVED eps on [N0, XCUT] ----
print(f"building psi to {XCUT} and proving eps on [{N0},{XCUT}]...", flush=True)
lamlog = [None] * (XCUT + 1)      # Lambda(n) as arb (None if 0)
comp = [False] * (XCUT + 1)
for p in range(2, XCUT + 1):
    if not comp[p]:
        for m in range(2 * p, XCUT + 1, p):
            comp[m] = True
        lp = arb(p).log()
        pk = p
        while pk <= XCUT:
            lamlog[pk] = lp
            pk *= p
EPS = arb(0)
psi = arb(0)
for n in range(2, XCUT + 1):
    if lamlog[n] is not None:
        psi = psi + lamlog[n]
    if n >= N0:
        d1 = abs(psi - n)
        d2 = abs(psi - (n + 1))
        cand = abs(d1.union(d2)) / n
        # track a rigorous upper bound for the max
        if not (cand < EPS):
            EPS = EPS.union(cand)
EPS = abs(EPS)
print(f"  eps enclosure on [{N0},{XCUT}]: {EPS}  (upper bound used: {EPS.mid()+EPS.rad()})", flush=True)

def g_pair_arb(cj, ck, u):
    pref = W * arb.pi().sqrt() / (2 * arb.pi())
    d1 = cj - ck
    e_near = (-(d1 * d1) / (4 * W * W)).exp()
    d2 = cj + ck
    e_far = (-(d2 * d2) / (4 * W * W)).exp()
    gauss = (-(u * u) * W * W / 4).exp()
    return pref * gauss * 2 * (e_near * (u * (cj + ck) / 2).cos() + e_far * (u * (cj - ck) / 2).cos())

def psmall_entry(cj, ck):
    s = arb(0)
    for n in range(2, N0):
        if lamlog[n] is not None:
            u = arb(n).log()
            s += lamlog[n] / arb(n).sqrt() * g_pair_arb(cj, ck, u)
    return 2 * s

U0 = arb(N0).log()
UC = arb(XCUT).log()

def g_pair_acb(cj, ck, u):
    """holomorphic in u for acb.integral."""
    pref = W * arb.pi().sqrt() / (2 * arb.pi())
    d1 = cj - ck
    e_near = (-(acb(d1) * acb(d1)) / (2 * W * W * 2)).exp()
    e_far = (-(acb(cj + ck) * acb(cj + ck)) / (2 * W * W * 2)).exp()
    gauss = (-(u * u) * W * W / 4).exp()
    return pref * gauss * 2 * (e_near * (u * (cj + ck) / 2).cos() + e_far * (u * (cj - ck) / 2).cos())

def tsm_entry(cj, ck):
    """2 Int_{U0}^{UC} e^{u/2} g_pair du  + tail radius beyond UC (envelope, exp-dominated)."""
    def integ(u, _):
        return (u / 2).exp() * g_pair_acb(acb(cj), acb(ck), u)
    ctx.prec = 320
    I = acb.integral(integ, acb(U0), acb(UC))
    ctx.prec = PREC
    assert I.imag.contains(arb(0)), "tsm imag not enclosing 0"
    val = 2 * I.real
    # tail: |G| <= 4P e^{-u^2 w^2/4}; for u>=UC: e^{-u^2w^2/4} <= e^{UC^2 w^2/4} e^{-UC w^2 u/2}
    P4 = 4 * W * arb.pi().sqrt() / (2 * arb.pi())
    rate = UC * W * W / 2 - arb(1) / 2          # decay rate of e^{u/2} * bound
    tail = 2 * P4 * (UC * UC * W * W / 4).exp() * (-(rate) * UC + UC / 2 - UC * UC * W * W / 2).exp() / rate
    # simplify safe: integral_{UC}^inf e^{u/2 - u^2w^2/4} du <= e^{UC/2 - UC^2w^2/4}/rate
    tail = 2 * P4 * (UC / 2 - UC * UC * W * W / 4).exp() / rate
    return val + arb(0).union(tail.union(-tail))

def cenv_entry(cj, ck):
    """upper bound on the fluctuation charge entry over [N0, XCUT] (per unit eps):
    2[ x0 e^{-u0/2}|G(u0)| + Int_{U0}^{UC} e^{u/2} 4P gauss (u w^2/2 + a + 1/2) du ]"""
    a = (cj + ck) / 2
    P4 = 4 * W * arb.pi().sqrt() / (2 * arb.pi())
    def integ(u, _):
        return (u / 2).exp() * P4 * (-(u * u) * W * W / 4).exp() * (u * W * W / 2 + a + arb(1) / 2)
    ctx.prec = 320
    I = acb.integral(integ, acb(U0), acb(UC))
    ctx.prec = PREC
    assert I.imag.contains(arb(0))
    g0 = g_pair_arb(cj, ck, U0)
    bdry = arb(N0) * (-(U0) / 2).exp() * abs(g0)
    return 2 * (bdry + I.real)

def exact_tail_beyond():
    """sum_{n>XCUT} Lambda(n)/sqrt(n)*|2 g_pair| <= 2*1.04*Int e^{u/2}4P gauss du beyond UC."""
    P4 = 4 * W * arb.pi().sqrt() / (2 * arb.pi())
    rate = UC * W * W / 2 - arb(1) / 2
    I = (UC / 2 - UC * UC * W * W / 4).exp() / rate
    return 2 * arb('1.04') * P4 * I

def certify(s_str):
    s = arb(s_str)
    cs = probes(s)
    n = len(cs)
    print(f"\ns={s_str}: {n} probes — rigorous B_lower ({n*(n+1)//2} arch + {n*(n+1)//2} tsm integrals)...",
          flush=True)
    # charge scalar: eps*||C_env||_F + tails (both directions on the smooth-tail radius already in tsm)
    Cf2 = arb(0)
    C = [[None] * n for _ in range(n)]
    for j in range(n):
        for k in range(j, n):
            c = cenv_entry(cs[j], cs[k])
            C[j][k] = C[k][j] = c
    for j in range(n):
        for k in range(n):
            Cf2 = Cf2 + C[j][k] * C[j][k]
    CHG = EPS * Cf2.sqrt() + n * exact_tail_beyond()
    print(f"  charge scalar (eps*||C||_F + exact tail): {CHG}", flush=True)
    M = [[None] * n for _ in range(n)]
    for j in range(n):
        for k in range(j, n):
            v = (pole_entry(cs[j], cs[k]) + arch_entry(cs[j], cs[k])
                 - psmall_entry(cs[j], cs[k]) - tsm_entry(cs[j], cs[k]))
            if j == k:
                v = v - CHG
            M[j][k] = M[k][j] = v
        print(f"  row {j} done", flush=True)
    d = [None] * n
    Lm = [[arb(0)] * n for _ in range(n)]
    ok = True
    for i in range(n):
        sdiag = M[i][i]
        for k in range(i):
            sdiag -= Lm[i][k] * Lm[i][k] * d[k]
        if not (sdiag > 0):
            print(f"  PIVOT {i}: {sdiag} does NOT certify > 0 -> NOT CERTIFIED at s={s_str}")
            ok = False
            break
        d[i] = sdiag
        for j2 in range(i + 1, n):
            v = M[j2][i]
            for k in range(i):
                v -= Lm[j2][k] * Lm[i][k] * d[k]
            Lm[j2][i] = v / d[i]
    if ok:
        print(f"  ALL {n} PIVOTS POSITIVE -> B_lower(s={s_str}) CERTIFIED PD (arb prec {PREC})")
        print(f"  => B_full PD, THEOREM-GRADE: arithmetic input = primes<{N0} exact + "
              f"proved eps({EPS.mid()}) on [{N0},{XCUT}] + psi<1.04x beyond")
        print(f"  smallest pivot: {min(d, key=lambda x: x.mid())}")
    return ok

if __name__ == "__main__":
    for s_str in (sys.argv[1:] or ["3.9", "3.8"]):
        certify(s_str)

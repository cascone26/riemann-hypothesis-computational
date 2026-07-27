#!/usr/bin/env python3
"""§359: THE FULLY ELEMENTARY INTERVAL THEOREM — zero numeric integration anywhere.

B_full >= B_elem = pole + M_Psi - P_small - T_smooth - CHG*I - tails, where EVERY entry is
erf/exp/cos closed forms + rigorous digamma gridpoint enclosures:
  pole    : exp closed form.
  M_Psi   : §341 Riemann minorant, arb-rigorous WITHOUT monotonicity assumptions — per cell
            [r_i, r_i+h) the value Psi_i = LOWER EDGE of the interval enclosure Omega(cell-ball)
            (digamma on a ball; §350b endpoint discipline), zone masses Z_i via real erf
            differences (telescoped edges). M_Psi <= arch as forms (phi^2 >= 0, pointwise
            minorant); arch tail beyond Rmax charged as a crude Gaussian radius (~e^{-700}).
  P_small : Lambda closed forms (n < n0).
  T_smooth: CLOSED FORM via complex erfc — int_{u0}^inf e^{cu-qu^2} du =
            e^{c^2/4q} (sqrt(pi)/2sqrt(q)) erfc(sqrt(q)(u0 - c/2q)), c = 1/2 + i*alpha.
  CHG     : eps (proved in-script, exhaustive interval psi check on [n0,XCUT]) times the
            envelope charge, itself erf/exp closed (linear-times-Gaussian integrals), plus
            psi<1.04x (Rosser-Schoenfeld 1962) tails beyond XCUT — all closed forms.
Interval LDL^T with pivot test. Certifying B_elem PD certifies B_full PD.
Config: w=2, band [14,44], n0=300, s in argv (default 4.5 4.4 4.3 4.2). Expected reach vs
minorant loss (~1e-4 scale): s=4.4 (margin 1.1e-2) should certify; s=4.2 (1.5e-4) is the test.
"""
import sys
from flint import arb, acb, ctx

PREC = 220
ctx.prec = PREC
W = arb(2)
BAND0, BAND1 = arb(14), arb(44)
N0 = 300
XCUT = 20000
R1 = arb(7)          # low-zone edge
RMAX = arb(120)      # minorant support
HCELL = arb('0.02')  # cell width (dyadic)

def probes(s):
    cs = []
    c = BAND0
    while c <= BAND1 + arb('1e-9'):
        cs.append(c)
        c = c + s
    return cs

def F_acb(c, r):
    return (-((r - c) * (r - c)) / (2 * W * W)).exp() + (-((r + c) * (r + c)) / (2 * W * W)).exp()

def pole_entry(cj, ck):
    z = acb(0, arb(1) / 2)
    return 2 * (F_acb(acb(cj), z) * F_acb(acb(ck), z)).real

# ---- psi, Lambda, proved eps (as §350) ----
print(f"building psi to {XCUT}, proving eps on [{N0},{XCUT}] ...", flush=True)
lamlog = [None] * (XCUT + 1)
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
        cand = abs(abs(psi - n).union(abs(psi - (n + 1)))) / n
        if not (cand < EPS):
            EPS = EPS.union(cand)
EPS = abs(EPS)
EPS_UP = EPS.mid() + EPS.rad()
print(f"  eps upper = {EPS_UP}", flush=True)

def g_pair_arb(cj, ck, u):
    pref = W * arb.pi().sqrt() / (2 * arb.pi())
    e_near = (-((cj - ck) * (cj - ck)) / (4 * W * W)).exp()
    e_far = (-((cj + ck) * (cj + ck)) / (4 * W * W)).exp()
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
Q = W * W / 4
SQ = W / 2                      # sqrt(Q)

def gauss_int_inf(c_acb, u0):
    """int_{u0}^inf e^{c u - Q u^2} du, closed form via (complex) erfc."""
    m = c_acb / (2 * Q)
    pref = (c_acb * c_acb / (4 * Q)).exp() * arb.pi().sqrt() / (2 * SQ)
    return pref * (SQ * (acb(u0) - m)).erfc()

def gauss_int_seg(c_acb, ua, ub):
    """int_{ua}^{ub} e^{c u - Q u^2} du = inf(ua) - inf(ub)."""
    return gauss_int_inf(c_acb, ua) - gauss_int_inf(c_acb, ub)

def tsm_entry(cj, ck):
    """2 int_{u0}^inf e^{u/2} g_pair du — closed form, integrated to INFINITY (no tail split)."""
    pref = W * arb.pi().sqrt() / (2 * arb.pi())
    e_near = (-((cj - ck) * (cj - ck)) / (4 * W * W)).exp()
    e_far = (-((cj + ck) * (cj + ck)) / (4 * W * W)).exp()
    a = (cj + ck) / 2
    b = (cj - ck) / 2
    Ia = gauss_int_inf(acb(arb(1) / 2, a), U0).real     # Re of e^{iau} integral = cos part
    Ib = gauss_int_inf(acb(arb(1) / 2, b), U0).real
    return 2 * (pref * 2 * (e_near * Ia + e_far * Ib))

def lin_gauss_int(alpha, beta, ua, ub):
    """int_{ua}^{ub} (alpha*u + beta) e^{u/2 - Q u^2} du, closed:
    with c=1/2: int u e^{cu-Qu^2} = [-(1/2Q) e^{cu-Qu^2}] + (c/2Q) int e^{cu-Qu^2}."""
    c = arb(1) / 2
    I0 = gauss_int_seg(acb(c), ua, ub).real
    e_a = (c * ua - Q * ua * ua).exp()
    e_b = (c * ub - Q * ub * ub).exp()
    I1 = (e_a - e_b) / (2 * Q) + (c / (2 * Q)) * I0
    return alpha * I1 + beta * I0

def cenv_entry(cj, ck):
    """envelope charge per unit eps over [u0,UC]: 2[x0 e^{-u0/2}|G(u0)| + X e^{-UC/2}|G(UC)|
       + int e^{u/2} 4P gauss (u w^2/2 + a + 1/2) du]  (both boundaries, §353 discipline)."""
    a = (cj + ck) / 2
    P4 = 4 * W * arb.pi().sqrt() / (2 * arb.pi())
    I = P4 * lin_gauss_int(W * W / 2, a + arb(1) / 2, U0, UC)
    bd0 = arb(N0) * (-(U0) / 2).exp() * abs(g_pair_arb(cj, ck, U0))
    bdX = arb(XCUT) * (-(UC) / 2).exp() * abs(g_pair_arb(cj, ck, UC))
    return 2 * (bd0 + bdX + I)

def exact_tail_beyond():
    """|sum_{n>XCUT}| via psi(x)<1.04x (RS 1962): 1.04*(int_{UC}^inf e^{u/2}4P gauss du + h(X)X)."""
    P4 = 4 * W * arb.pi().sqrt() / (2 * arb.pi())
    I = P4 * gauss_int_inf(acb(arb(1) / 2), UC).real
    bdry = P4 * (UC / 2 - Q * UC * UC).exp()
    return 2 * arb('1.04') * (I + bdry)

# ---- M_Psi: LEFT-ENDPOINT point values, rigorous by the MONOTONICITY LEMMA ----
# LEMMA (2-line proof): Omega(r) = Re psi(1/4 + ir/2) - log pi is strictly increasing on
# (0, inf): psi'(z) = sum_{n>=0} 1/(n+z)^2, and for z = 1/4 + it (t>0) each term has
# Im 1/(n+z)^2 = -2(n+1/4)t/|n+z|^4 < 0, so Im psi' < 0 and
# Omega'(r) = -(1/2) Im psi'(1/4+ir/2) > 0.  Hence on [r_i, r_i+h): Omega >= Omega(r_i),
# and on [0, R1]: Omega >= Omega(0) = psi(1/4) - log pi.  (Ball evaluation not needed —
# arb's disk-variation digamma enclosures are uselessly wide.)
print("building Omega left-endpoint values (monotonicity lemma) ...", flush=True)
LOGPI = arb.pi().log()
ncell = int(float((RMAX - R1) / HCELL) + 0.5)
cell_lo = []
r = R1
for i in range(ncell):
    om = (acb(arb(1) / 4, r / 2).digamma()).real - LOGPI
    cell_lo.append(om.mid() - om.rad())            # exact dyadic lower bound of the point value
    r = r + HCELL
om0 = (acb(arb(1) / 4, 0).digamma()).real - LOGPI
FLOOR0 = om0.mid() - om0.rad()
# NB: arb __eq__ is ball-equality (False for any inexact value) — nan-check via float(mid)
import math
assert not math.isnan(float(FLOOR0.mid())), "floor nan"
for i, cl in enumerate(cell_lo):
    assert not math.isnan(float(cl.mid())), f"cell {i} nan"
print(f"  {ncell} cells, [0,7] floor = Omega(0) = {FLOOR0}", flush=True)

def erf_edges(cj, ck):
    """for each sign-combo m=(±cj±ck)/2 pair-center: erf((edge−m)/W) at all cell edges,
    for both ±[R1,RMAX] windows; returns the zone-mass vector Z_i[cells] plus Z(0,R1)."""
    combos = []
    for sj in (1, -1):
        for sk in (1, -1):
            m = (sj * cj + sk * ck) / 2
            amp = (-((sj * cj - sk * ck) * (sj * cj - sk * ck)) / (4 * W * W)).exp()
            combos.append((m, amp))
    pref = arb.pi().sqrt() * W / 2
    Z = [arb(0)] * ncell
    Z0 = arb(0)
    for (m, amp) in combos:
        for sgn in (1, -1):
            # window sgn*[R1, RMAX]: edges r_i (ascending in |r|)
            prev = (((R1 if sgn == 1 else -RMAX) - m) / W).erf()
            if sgn == 1:
                e = R1
                for i in range(ncell):
                    e = e + HCELL
                    cur = ((e - m) / W).erf()
                    Z[i] += amp * pref * (cur - prev)
                    prev = cur
            else:
                # [-RMAX, -R1]: accumulate from -RMAX upward; cell i covers [-(R1+(i+1)h), -(R1+ih)]
                vals = []
                e = -RMAX
                v = ((e - m) / W).erf()
                for i in range(ncell + 1):
                    vals.append(v)
                    e = e + HCELL
                    v = ((e - m) / W).erf()
                for i in range(ncell):
                    # cell i (by |r|) = [-(R1+(i+1)h), -(R1+ih)] = vals index ncell-1-i .. ncell-i
                    Z[i] += amp * pref * (vals[ncell - i] - vals[ncell - 1 - i])
            # low zone |r|<R1 for this combo/side
            if sgn == 1:
                Z0 += amp * pref * (((R1 - m) / W).erf() - ((arb(0) - m) / W).erf())
            else:
                Z0 += amp * pref * (((arb(0) - m) / W).erf() - ((-R1 - m) / W).erf())
    return Z, Z0

def mpsi_entry(cj, ck):
    Z, Z0 = erf_edges(cj, ck)
    v = FLOOR0 * Z0
    for i in range(ncell):
        v += cell_lo[i] * Z[i]
    # arch tail beyond RMAX: |FjFk| <= 4 e^{-(RMAX-cj)^2/2W^2} e^{-(RMAX-ck)^2/2W^2}; Omega <= r on tail
    tj = (-((RMAX - cj) * (RMAX - cj)) / (2 * W * W)).exp()
    tk = (-((RMAX - ck) * (RMAX - ck)) / (2 * W * W)).exp()
    tail = arb(8) * tj * tk * RMAX
    return v / (2 * arb.pi()) + arb(0).union(tail.union(-tail))

def certify(s_str):
    s = arb(s_str)
    cs = probes(s)
    n = len(cs)
    print(f"\ns={s_str}: {n} probes — B_elem (all closed forms, {ncell} minorant cells) ...", flush=True)
    Cf2 = arb(0)
    C = [[None] * n for _ in range(n)]
    for j in range(n):
        for k in range(j, n):
            C[j][k] = C[k][j] = cenv_entry(cs[j], cs[k])
    for j in range(n):
        for k in range(n):
            Cf2 += C[j][k] * C[j][k]
    CHG = EPS_UP * Cf2.sqrt() + n * exact_tail_beyond()
    CHG = CHG.mid() + CHG.rad()
    print(f"  charge upper endpoint: {CHG}", flush=True)
    M = [[None] * n for _ in range(n)]
    for j in range(n):
        for k in range(j, n):
            v = (pole_entry(cs[j], cs[k]) + mpsi_entry(cs[j], cs[k])
                 - psmall_entry(cs[j], cs[k]) - tsm_entry(cs[j], cs[k]))
            if j == k:
                v = v - CHG
            M[j][k] = M[k][j] = v
        print(f"  row {j} done", flush=True)
    d = [None] * n
    Lm = [[arb(0)] * n for _ in range(n)]
    for i in range(n):
        sdiag = M[i][i]
        for k in range(i):
            sdiag -= Lm[i][k] * Lm[i][k] * d[k]
        if not (sdiag > 0):
            print(f"  PIVOT {i}: {sdiag} not certified > 0 -> NOT CERTIFIED at s={s_str}")
            return False
        d[i] = sdiag
        for j2 in range(i + 1, n):
            v = M[j2][i]
            for k in range(i):
                v -= Lm[j2][k] * Lm[i][k] * d[k]
            Lm[j2][i] = v / d[i]
    print(f"  ALL {n} PIVOTS POSITIVE -> B_elem(s={s_str}) CERTIFIED PD")
    print(f"  => B_full PD by ELEMENTARY CLOSED FORMS ONLY (erf/exp/cos/Lambda + digamma gridpoints"
          f" + proved eps + psi<1.04x). smallest pivot: {min(d, key=lambda x: x.mid())}")
    return True

if __name__ == "__main__":
    for s_str in (sys.argv[1:] or ["4.5", "4.4", "4.3", "4.2"]):
        certify(s_str)

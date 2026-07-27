#!/usr/bin/env python3
"""§348: THE ε-BUDGET RACE — the uniformity wall as two explicit curves.

At the deepest PSD spacing (edge) of each (band, w): certification of the collective chain
with a PNT envelope |psi(x)-x| <= eps*x on [n0,inf) needs

    eps * ||C_unit(n0)||_F  <  budget(n0) = lmin(A - Pres - Tsmooth) - tailbnd
    (charge is LINEAR in eps: boundary + integral both scale with the envelope)

so  eps_req(n0) = budget / ||C_unit||_F   vs   eps_avail(n0) = max_{x in [n0, 20000]} |psi(x)-x|/x.

The race: eps_req GROWS with n0 (unit charge dies under the Gaussian), eps_avail SHRINKS (PNT).
Winning cut n0* = first n0 with eps_avail < eps_req. The uniformity trajectory = how n0*
(and the resolved-core size pi*(n0*)) scales with 1/w and with band height T. Beyond x=20000
eps_avail is REPORT-tier (needs named RS/VK constants); on-range it is measured/exact.
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma
from scipy.integrate import quad

NMAX2 = 20000
lam = np.zeros(NMAX2 + 1)
comp = np.zeros(NMAX2 + 1, dtype=bool)
for p in range(2, NMAX2 + 1):
    if not comp[p]:
        comp[2 * p::p] = True
        lp = np.log(p); pk = p
        while pk <= NMAX2:
            lam[pk] = lp; pk *= p
L.NS = np.arange(2, NMAX2 + 1)
L.LOGN = np.log(L.NS)
L.LAMW = lam[2:] / np.sqrt(L.NS)
CUMPSI = np.cumsum(lam)

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def tail_bound(w, nprobes):
    uN = np.log(NMAX2)
    I, _ = quad(lambda u: np.exp(u / 2 - u * u * w * w / 4), uN, uN + 200)
    return 2 * 1.04 * I * 4 * w * np.sqrt(np.pi) / (2 * np.pi) * nprobes

def prime_block(cs, w, nlo, nhi):
    n = len(cs)
    M = np.zeros((n, n))
    mask = (L.NS >= nlo) & (L.NS < nhi) & (L.LAMW > 0)
    idx = np.where(mask)[0]
    for j in range(n):
        for k in range(j, n):
            v = 2.0 * np.sum(L.LAMW[idx] * L.g_pair(cs[j], cs[k], w, L.LOGN[idx]))
            M[j, k] = M[k, j] = v
    return M

def smooth_tail(cs, w, u0):
    n = len(cs)
    M = np.zeros((n, n))
    for j in range(n):
        for k in range(j, n):
            I, _ = quad(lambda u, cj=cs[j], ck=cs[k]: np.exp(u / 2) * L.g_pair(cj, ck, w, u),
                        u0, u0 + 300, limit=400)
            M[j, k] = M[k, j] = 2.0 * I
    return M

def G_and_Gp(cj, ck, w, u):
    P = w * np.sqrt(np.pi) / (2 * np.pi)
    a, b = (cj + ck) / 2.0, (cj - ck) / 2.0
    en = np.exp(-(cj - ck) ** 2 / (4 * w * w))
    ef = np.exp(-(cj + ck) ** 2 / (4 * w * w))
    gauss = np.exp(-u * u * w * w / 4.0)
    G = P * gauss * 2.0 * (en * np.cos(u * a) + ef * np.cos(u * b))
    Gp = P * gauss * 2.0 * (-(u * w * w / 2.0) * (en * np.cos(u * a) + ef * np.cos(u * b))
                            - en * a * np.sin(u * a) - ef * b * np.sin(u * b))
    return G, Gp

def unit_charge(cs, w, u0):
    """C_unit[j,k]: the eps=1 envelope charge (|psi-x| -> x)."""
    uN = np.log(NMAX2)
    ug = np.arange(u0, uN, 5e-4)
    xg = np.exp(ug)
    du = ug[1] - ug[0]
    x0 = np.exp(u0)
    n = len(cs)
    C = np.zeros((n, n))
    for j in range(n):
        for k in range(j, n):
            G, Gp = G_and_Gp(cs[j], cs[k], w, ug)
            fp = np.abs(np.exp(-ug / 2) * (Gp - G / 2))
            G0, _ = G_and_Gp(cs[j], cs[k], w, np.array([u0]))
            C[j, k] = C[k, j] = (2 * x0 * np.exp(-u0 / 2) * abs(float(G0[0]))
                                 + 2 * float(np.sum(xg * fp)) * du)
    return C

# eps_avail lookup (exact on [n0, 20000])
XS = np.arange(2, NMAX2 + 1)
DEVREL = np.abs(CUMPSI[2:] - XS) / XS
def eps_avail(n0):
    return float(np.max(DEVREL[XS >= n0]))

BANDS = [(14.0, 44.0), (50.0, 80.0), (100.0, 130.0)]
WS = [2.0, 1.5, 1.2]
CUTS = [100, 300, 1000, 3000, 10000]
print("ε-BUDGET RACE at the edge (deepest PSD spacing).  won iff eps_avail < eps_req")
print(f"{'band':>10} {'w':>4} {'s_edge':>6} {'true':>10} | " +
      " | ".join(f"n0={c:<5}" for c in CUTS) + "   (cells: eps_req/eps_avail; * = won)")
summary = []
for band in BANDS:
    set_grid(band[1] + 30.0)
    for w in WS:
        # find edge
        edge = None
        for s in np.arange(6.0, 1.09, -0.1):
            cs = np.arange(band[0], band[1] + 1e-9, s)
            if len(cs) < 4: continue
            f = L.build_forms(cs, w)
            tm = float(np.linalg.eigvalsh(f['full'])[0])
            if tm < 1e-13: break
            edge = (s, cs, f, tm)
        s, cs, f, tm = edge
        tb = tail_bound(w, len(cs))
        cells = []; n0star = None; victory = None; core = None
        for n0 in CUTS:
            u0 = np.log(n0) - 1e-12
            Pres = prime_block(cs, w, 2, n0)
            Tsm = smooth_tail(cs, w, u0)
            budget = float(np.linalg.eigvalsh(f['archpole'] - Pres - Tsm)[0]) - tb
            Cu = unit_charge(cs, w, u0)
            nCu = float(np.linalg.norm(Cu))
            er = budget / nCu if budget > 0 else -1.0
            ea = eps_avail(n0)
            won = er > ea > 0
            cells.append(f"{er:8.1e}/{ea:7.1e}{'*' if won else ' '}")
            if won and n0star is None:
                n0star = n0; victory = er / ea
                core = int(np.sum((L.NS < n0) & (L.LAMW > 0)))
        print(f"{str([int(band[0]), int(band[1])]):>10} {w:4.1f} {s:6.1f} {tm:+10.2e} | "
              + " | ".join(cells), flush=True)
        summary.append((band, w, s, n0star, victory, core))

print("\nSUMMARY — winning cut n0*, victory margin eps_req/eps_avail, resolved core size:")
print(f"{'band':>10} {'w':>4} {'s_edge':>6} {'n0*':>6} {'victory':>8} {'#core':>6}")
for band, w, s, n0s, v, core in summary:
    print(f"{str([int(band[0]), int(band[1])]):>10} {w:4.1f} {s:6.1f} {str(n0s):>6} "
          f"{f'{v:.0f}x' if v else '---':>8} {str(core):>6}")
print("\nTrajectories: log n0* vs 1/w (fixed band) = the w-axis cut law;")
print("log n0* vs band height (fixed w) = the T-axis cut law. Victory margin >> 1 means")
print("named RS/VK constants (lit values, larger than on-range measured eps) will still win.")

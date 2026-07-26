#!/usr/bin/env python3
"""§347: A-PRIORI PRICING of the fluctuation charge ||DT|| (the §346 named ingredient).

DT[j,k] = 2 Int_{x0}^inf x^{-1/2} G(log x) d(psi(x) - x),  G = g_pair(cj,ck,w,u).
Partial summation (f(u) = e^{-u/2} G(u), x = e^u):
  |DT[j,k]| <= 2|psi(x0)-x0| e^{-u0/2}|G(u0)|
             + 2 Int_{u0}^{uN} |psi(e^u)-e^u| |f'(u)| du        (exact psi to NMAX)
             + tail beyond NMAX (psi envelope, Gaussian-killed)
  f'(u) = e^{-u/2} (G'(u) - G(u)/2),  G' analytic.

Two envelopes per entry:
  EXACT : |psi(x)-x| pointwise (the honest partial-summation cost; the gap to the measured
          ||DT|| isolates the loss from taking |.| — i.e. the cancellation INSIDE the
          fluctuation that pointwise partial summation cannot see)
  CONST : |psi(x)-x| <= eps*x with eps = max over [x0, NMAX] measured (stand-in for an
          explicit RS-type constant — same shape as the fully a-priori version; the
          named-constant substitution is bookkeeping, not structure)

Charge matrix C >= |DT| entrywise  =>  ||DT||_2 <= ||C||_F.  Chain:
  cert = lmin(A - Pres - Tsmooth) - ||C||_F - tailbnd.
RACE: charge loss (orders) vs self-averaging gain (~3 orders per cut step) — push n0 up to
5000 and see if the a-priori chain still reaches mid-plunge with a few hundred resolved primes.
"""
import numpy as np
import l2_edge_sweep as L
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

def charge_matrices(cs, w, u0):
    """returns (C_exact, C_const, eps_used): entrywise partial-summation bounds on |DT|."""
    uN = np.log(NMAX2)
    ug = np.arange(u0, uN, 5e-4)
    xg = np.exp(ug)
    psig = CUMPSI[np.floor(xg).astype(int)]
    dev = np.abs(psig - xg)                       # exact |psi - x| on the grid
    eps = float(np.max(dev / xg))
    x0 = np.exp(u0)
    dev0 = abs(float(CUMPSI[int(np.floor(x0))]) - x0)
    n = len(cs)
    Ce = np.zeros((n, n)); Cc = np.zeros((n, n))
    du = ug[1] - ug[0]
    for j in range(n):
        for k in range(j, n):
            G, Gp = G_and_Gp(cs[j], cs[k], w, ug)
            fp = np.abs(np.exp(-ug / 2) * (Gp - G / 2))
            G0, _ = G_and_Gp(cs[j], cs[k], w, np.array([u0]))
            bdry = 2 * dev0 * np.exp(-u0 / 2) * abs(float(G0[0]))
            Ce[j, k] = Ce[k, j] = bdry + 2 * float(np.sum(dev * fp)) * du
            Cc[j, k] = Cc[k, j] = bdry + 2 * float(np.sum(eps * xg * fp)) * du
    return Ce, Cc, eps

BAND = (14.0, 44.0)
print("A-PRIORI PRICING of ||DT|| + the cut race.  cert = lmin(M) - charge - tailbnd")
print(f"{'w':>4} {'n0':>5} {'#res':>5} {'meas||DT||':>10} {'||Ce||F':>9} {'loss_e':>7} "
      f"{'||Cc||F':>9} {'eps':>6} {'reach_meas':>10} {'reach_exact':>11} {'reach_const':>11}")
for w in (2.0, 1.5, 1.2):
    tb = None
    for n0 in (100, 500, 2000, 5000):
        u0 = np.log(n0) - 1e-12
        nres = int(np.sum((L.NS < n0) & (L.LAMW > 0)))
        # reference s for the measured ||DT|| (mid-plunge probe spacing)
        reach = {'m': None, 'e': None, 'c': None}
        printed = False
        for s in np.arange(6.0, 1.19, -0.2):
            cs = np.arange(BAND[0], BAND[1] + 1e-9, s)
            if len(cs) < 4: continue
            f = L.build_forms(cs, w)
            true_m = float(np.linalg.eigvalsh(f['full'])[0])
            if true_m < 1e-13: break
            if tb is None: tb = tail_bound(w, len(cs))
            Pres = prime_block(cs, w, 2, n0)
            Tsm = smooth_tail(cs, w, u0)
            DT = (f['prime'] - Pres) - Tsm
            nDT = float(np.linalg.norm(DT, 2))
            Ce, Cc, eps = charge_matrices(cs, w, u0)
            nCe, nCc = float(np.linalg.norm(Ce)), float(np.linalg.norm(Cc))
            lmin = float(np.linalg.eigvalsh(f['archpole'] - Pres - Tsm)[0])
            if lmin - nDT - tb > 0: reach['m'] = s
            if lmin - nCe - tb > 0: reach['e'] = s
            if lmin - nCc - tb > 0: reach['c'] = s
            if not printed and s < 4.7:   # record the mid-plunge row for the table
                row = (nDT, nCe, nCc, eps)
                printed = True
        nDT, nCe, nCc, eps = row
        loss = nCe / nDT if nDT > 0 else float('inf')
        print(f"{w:4.1f} {n0:5d} {nres:5d} {nDT:10.2e} {nCe:9.2e} {loss:7.0f} "
              f"{nCc:9.2e} {eps:6.3f} {str(reach['m']):>10} {str(reach['e']):>11} "
              f"{str(reach['c']):>11}", flush=True)
print("\nReading: loss_e = pure |.|-loss (cancellation inside the fluctuation that pointwise")
print("partial summation cannot see). If reach_exact/const still hits ~4.4-4.6 at n0=2000-5000,")
print("the a-priori collective chain closes mid-plunge with a few hundred resolved primes;")
print("if reach dies at every cut, the charge loss beats self-averaging and mean-square")
print("(large-sieve) pricing replaces pointwise partial summation as the named ingredient.")

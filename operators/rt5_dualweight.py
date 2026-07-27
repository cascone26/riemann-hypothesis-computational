#!/usr/bin/env python3
"""§357: the PRIME<->POLE dual-weight probe (the §354/§356 named next face).

At each band, descending s toward the edge (w=2): v = bottom eigenvector of A = arch+pole;
measure the minimizer projections p_pole = <v, pole v>, p_arch = <v, arch v>,
p_n* = |<v, P_n* v>| (resonant), and sum_res = sum over projection-eff10 of |<v,P_n v>|.
Questions (stated before the run): (i) does p_n*/p_pole converge to a clean band-independent
law approaching the edge? (ii) does the TOTAL major-arc projection track the pole projection
(the involution-pairing made quantitative)? A stable ratio => the resonant term admits a dual
representation in the exactly-known pole term. No pre-registered numeric (exploratory probe —
labeled as such); whatever law appears must then be derived before it charges anything.
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

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

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def prime_terms_w(cs, w, thresh=1e-16):
    out = {}
    for i, n in enumerate(L.NS):
        if L.LAMW[i] == 0.0: continue
        u = L.LOGN[i]
        if L.LAMW[i] * np.exp(-u * u * w * w / 4) < thresh: continue
        M = np.zeros((len(cs), len(cs)))
        for j in range(len(cs)):
            for k in range(j, len(cs)):
                v = 2.0 * L.LAMW[i] * L.g_pair(cs[j], cs[k], w, u)
                M[j, k] = M[k, j] = v
        out[int(n)] = M
    return out

w = 2.0
BANDS = [(14.0, 44.0), (50.0, 80.0), (100.0, 130.0), (200.0, 230.0)]
print("PRIME<->POLE dual weight, w=2 (exploratory probe — no charge until derived)")
print(f"{'band':>10} {'s':>5} {'true':>10} {'p_pole':>10} {'p_arch':>9} {'n*':>4} {'p_n*':>10} "
      f"{'p_n*/p_arch':>11} {'sum_res':>10} {'sumres/arch':>11}")
for band in BANDS:
    set_grid(band[1] + 30.0)
    rows = 0
    for s in np.arange(6.0, 1.09, -0.1):
        cs = np.arange(band[0], band[1] + 1e-9, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        tm = float(np.linalg.eigvalsh(f['full'])[0])
        if tm < 1e-13: break
        laste = (s, cs, f, tm)
    # walk the last 4 spacings above the edge
    s_edge = laste[0]
    for s in (s_edge + 0.3, s_edge + 0.2, s_edge + 0.1, s_edge):
        cs = np.arange(band[0], band[1] + 1e-9, s)
        f = L.build_forms(cs, w)
        tm = float(np.linalg.eigvalsh(f['full'])[0])
        # pole matrix
        z = 0.5j
        fp = np.exp(-((z - cs) ** 2) / (2 * w * w)) + np.exp(-((z + cs) ** 2) / (2 * w * w))
        pole = 2.0 * np.real(np.outer(fp, fp))
        arch = f['archpole'] - pole
        lA, UA = np.linalg.eigh(f['archpole'])
        v = UA[:, 0]
        p_pole = float(v @ pole @ v)
        p_arch = float(v @ arch @ v)
        terms = prime_terms_w(cs, w)
        proj = {n: abs(float(v @ terms[n] @ v)) for n in terms}
        order = sorted(proj, key=lambda n: -proj[n])
        nstar = order[0]
        top = proj[nstar]
        sumres = sum(proj[n] for n in order if proj[n] >= top / 10)
        print(f"{str([int(band[0]), int(band[1])]):>10} {s:5.1f} {tm:+10.2e} {p_pole:10.3e} "
              f"{p_arch:9.3e} {nstar:4d} {top:10.3e} {top/p_arch:11.3e} {sumres:10.3e} "
              f"{sumres/p_arch:11.3e}", flush=True)

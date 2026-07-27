#!/usr/bin/env python3
"""§354b: Q4 calibration REDONE with the RIGHT resonance object.

§354a caught a selection-notion conflation: raw term norm ||P_n|| is maximized by n=2 always
(Lambda/sqrt(n) weight); the §337 Bragg law log n* ~ 2pi/s lives in the MINIMIZER PROJECTION
p_n = |<v, P_n v>|, v = bottom eigenvector of A = arch+pole. Redo: per band at the w=2 edge,
print the projection spectrum's top terms, n*_proj, and run the involution-closure test on the
projection-eff10 set. Also print projections at criticality vs the three candidate laws.
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

BANDS = [(14.0, 44.0), (50.0, 80.0), (100.0, 130.0), (200.0, 230.0)]
w = 2.0
print("Q4 REDO — resonance by MINIMIZER PROJECTION p_n = |<v, P_n v>|, v = argmin eig(A)")
print(f"{'band':>10} {'s_edge':>6} {'n*proj':>6} {'log n*':>7} {'2pi/s':>6} {'log(gc/2pi)':>11} "
      f"{'top-5 by projection':>40}")
results = []
for band in BANDS:
    set_grid(band[1] + 30.0)
    edge = None
    for s in np.arange(6.0, 1.09, -0.1):
        cs = np.arange(band[0], band[1] + 1e-9, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        tm = float(np.linalg.eigvalsh(f['full'])[0])
        if tm < 1e-13: break
        edge = (s, cs, f, tm)
    s, cs, f, tm = edge
    terms = prime_terms_w(cs, w)
    lA, UA = np.linalg.eigh(f['archpole'])
    v = UA[:, 0]
    proj = {n: abs(float(v @ terms[n] @ v)) for n in terms}
    order = sorted(proj, key=lambda n: -proj[n])
    nstar = order[0]
    gc = (band[0] + band[1]) / 2
    top5 = ", ".join(f"{n}:{proj[n]:.2e}" for n in order[:5])
    print(f"{str([int(band[0]), int(band[1])]):>10} {s:6.1f} {nstar:6d} {np.log(nstar):7.3f} "
          f"{2*np.pi/s:6.3f} {np.log(gc/(2*np.pi)):11.3f} {top5:>40}", flush=True)
    results.append((band, s, gc, proj, order, terms))

print("\nInvolution closure on projection-eff10 (images under n -> gc/(2pi n)):")
for band, s, gc, proj, order, terms in results:
    top = proj[order[0]]
    eff10 = [n for n in order if proj[n] >= top / 10]
    eff100set = set(n for n in order if proj[n] >= top / 100)
    rows = []
    for n in eff10[:6]:
        img = gc / (2 * np.pi * n)
        if img < 1.5:
            rows.append(f"{n}->pole-side"); continue
        cand = min(terms.keys(), key=lambda m: abs(np.log(m) - np.log(img)))
        du = abs(np.log(cand) - np.log(img))
        inside = cand in eff100set and du < np.pi / (2 * s)
        rows.append(f"{n}->{img:.1f}~{cand} {'IN' if inside else 'out'}")
    print(f"  {str([int(band[0]), int(band[1])]):>10} eff10={eff10[:6]}: " + " | ".join(rows))

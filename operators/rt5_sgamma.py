#!/usr/bin/env python3
"""§368: the s<->gamma calibration (grounding Q4, made quantitative).

Measured Bragg law: log n* ~ 2pi/s (projection object, §354). At the Nyquist/edge spacing
s_edge(T) ~ 1/rho(T) = 2pi/log(T/2pi), this predicts log n*_edge ~ 2pi/s_edge ~ log(T/2pi),
i.e. n*_edge ~ T/2pi. Test directly: at each band's edge, is n*_proj ~ T/(2pi)? If yes, the
critical resonant scale IS the RS scale gamma/2pi after all — NOT as a code artifact (§343
retraction stands: the code never computed it) but as a DERIVED consequence of Bragg + Nyquist,
now measured. This resolves Q4's calibration and re-tests the involution n -> (T/2pi)/n on the
correct (projection) resonant set.
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

NMAX2 = 20000
lam = np.zeros(NMAX2 + 1); comp = np.zeros(NMAX2 + 1, dtype=bool)
for p in range(2, NMAX2 + 1):
    if not comp[p]:
        comp[2 * p::p] = True
        lp = np.log(p); pk = p
        while pk <= NMAX2:
            lam[pk] = lp; pk *= p
L.NS = np.arange(2, NMAX2 + 1); L.LOGN = np.log(L.NS); L.LAMW = lam[2:] / np.sqrt(L.NS)

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
                M[j, k] = M[k, j] = 2.0 * L.LAMW[i] * L.g_pair(cs[j], cs[k], w, u)
        out[int(n)] = M
    return out

w = 2.0
CENTERS = [29, 45, 65, 90, 115, 150, 200, 260, 340, 440]
WID = 30.0
print("s<->gamma: is the projection resonant scale n*_edge ~ T/(2pi)?  (Bragg + Nyquist prediction)")
print(f"{'Tc':>5} {'s_edge':>6} {'2pi/s':>6} {'log(T/2pi)':>10} {'n*_proj':>7} {'T/2pi':>6} {'n*/(T/2pi)':>10}")
rows = []
for Tc in CENTERS:
    b0, b1 = Tc - WID / 2, Tc + WID / 2
    set_grid(b1 + 30.0)
    edge = None
    for s in np.arange(6.0, 1.05, -0.05):
        cs = np.arange(b0, b1 + 1e-9, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        tm = float(np.linalg.eigvalsh(f['full'])[0])
        if tm < 1e-13: break
        edge = (s, cs, f)
    s, cs, f = edge
    terms = prime_terms_w(cs, w)
    lA, UA = np.linalg.eigh(f['archpole'])
    v = UA[:, 0]
    proj = {n: abs(float(v @ terms[n] @ v)) for n in terms}
    nstar = max(proj, key=lambda n: proj[n])
    T2pi = Tc / (2 * np.pi)
    rows.append((Tc, s, nstar, T2pi))
    print(f"{Tc:5d} {s:6.2f} {2*np.pi/s:6.3f} {np.log(Tc/(2*np.pi)):10.3f} {nstar:7d} {T2pi:6.1f} "
          f"{nstar/T2pi:10.3f}", flush=True)

ratios = [r[2] / r[3] for r in rows]
print(f"\nn*/(T/2pi): mean {np.mean(ratios):.3f}, sd {np.std(ratios):.3f}")
print("verdict: if mean~O(1) and sd small, the RS scale gamma/2pi IS the critical resonant scale")
print("(derived from Bragg+Nyquist, not assumed) — Q4 calibration pinned; §343 retraction refined:")
print("the ANNOTATION was wrong to state as measured, but the RELATION is true and now IS measured.")
# involution self-consistency: does n* pair with T/2pi / n* ~ 1 (pole side)?
print(f"\ninvolution n -> (T/2pi)/n* on the resonant scale:")
for Tc, s, nstar, T2pi in rows:
    print(f"  T={Tc}: n*={nstar} -> {T2pi/nstar:.2f}  ({'pole/leading side' if T2pi/nstar < 3 else 'prime side'})")

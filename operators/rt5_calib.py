#!/usr/bin/env python3
"""§354: Q4 CALIBRATION + Q2 FACTORIZATION.

Q4 (s<->gamma at criticality, settles Dirac's involution branch): at each band's edge,
measure n* (top-norm term) and compare log n* against: 2pi/s_edge (Bragg), log(gamma_c/2pi)
(RS-endpoint reading), (1/2)log(gamma_c/2pi) (RS-crossover/fixed-point reading). Then test
closure of the eff10 major-arc set under the involution n -> gamma_c/(2pi n): images inside
the eff100 set => closed (pairing gain); images off => sieve absorbs. Either branch banks.

Q2 (derive R ~ 40/w^2): factor R = (#occupied Bragg slots) x (couplings per slot) at fixed
margin depth 1e-8 across (band, w). Hypothesis registered BEFORE the run: slots ~ 1/w and
couplings-per-slot ~ 1/w carry the quadratic; the T-growth of R (15->19->24 at w=2, §344)
should sit in ONE factor (predict: couplings-per-slot via the resonant scale e^{2pi/s}
climbing into denser primes as s_edge shrinks with T).
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

def gap_parts(M, Q, k=5):
    n = M.shape[0]; k = min(k, n - 1)
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :k], UM[:, k:]
    top = float(np.linalg.eigvalsh(np.diag(lM[:k]) - Ub.T @ Q @ Ub)[0])
    q = float(np.linalg.norm(Ub.T @ Q @ Ur, 2))
    pr = float(np.linalg.eigvalsh(Ur.T @ Q @ Ur)[-1])
    return float(np.linalg.eigvalsh(np.array([[top, -q], [-q, lM[k] - pr]]))[0])

BANDS = [(14.0, 44.0), (50.0, 80.0), (100.0, 130.0), (200.0, 230.0)]

print("=" * 100)
print("Q4 — s<->gamma calibration at the edge (w=2)")
print("=" * 100)
w = 2.0
print(f"{'band':>10} {'s_edge':>6} {'n*':>5} {'log n*':>7} {'2pi/s':>6} {'log(gc/2pi)':>11} "
      f"{'0.5log(gc/2pi)':>14}")
calib = []
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
    norms = {n: float(np.linalg.norm(terms[n], 2)) for n in terms}
    order = sorted(norms, key=lambda n: -norms[n])
    nstar = order[0]
    gc = (band[0] + band[1]) / 2
    print(f"{str([int(band[0]), int(band[1])]):>10} {s:6.1f} {nstar:5d} {np.log(nstar):7.3f} "
          f"{2*np.pi/s:6.3f} {np.log(gc/(2*np.pi)):11.3f} {0.5*np.log(gc/(2*np.pi)):14.3f}")
    calib.append((band, s, cs, f, terms, norms, order, gc))

print("\nInvolution closure test: images of eff10 set under n -> gc/(2pi n)")
for band, s, cs, f, terms, norms, order, gc in calib:
    top = norms[order[0]]
    eff10 = [n for n in order if norms[n] >= top / 10]
    eff100set = set(n for n in order if norms[n] >= top / 100)
    rows = []
    for n in eff10:
        img = gc / (2 * np.pi * n)
        # nearest prime power in log-distance
        if img < 1.5:
            rows.append(f"{n}->img {img:.2f} (<2: maps to POLE/leading side)")
            continue
        cand = min(terms.keys(), key=lambda m: abs(np.log(m) - np.log(img)))
        du = abs(np.log(cand) - np.log(img))
        inside = cand in eff100set and du < np.pi / (2 * s)
        rows.append(f"{n}->img {img:.2f}~pp {cand} (du={du:.2f}) {'IN eff100' if inside else 'minor arc'}")
    print(f"  {str([int(band[0]), int(band[1])]):>10}: " + " | ".join(rows))

print("\n" + "=" * 100)
print("Q2 — factorization R = slots x couplings-per-slot at margin depth 1e-8")
print("=" * 100)
print(f"{'band':>10} {'w':>4} {'s@1e-8':>7} {'R(k5)':>6} {'slots':>6} {'cps':>5} "
      f"{'slots*cps':>9} {'Rw2':>6} {'eff100':>6}")
for band in BANDS[:3]:
    set_grid(band[1] + 30.0)
    for w2 in (2.5, 2.0, 1.5):
        best = None
        for s in np.arange(6.0, 1.09, -0.05):
            cs = np.arange(band[0], band[1] + 1e-9, s)
            if len(cs) < 4: continue
            f = L.build_forms(cs, w2)
            tm = float(np.linalg.eigvalsh(f['full'])[0])
            if tm < 1e-13: break
            d = abs(np.log10(max(tm, 1e-300)) + 8)
            if best is None or d < best[0]: best = (d, s, cs, f, tm)
        _, s, cs, f, tm = best
        terms = prime_terms_w(cs, w2)
        norms = {n: float(np.linalg.norm(terms[n], 2)) for n in terms}
        order = sorted(norms, key=lambda n: -norms[n])
        A, P = f['archpole'], f['prime']
        Pres = np.zeros_like(A); R = None
        for cnt, n in enumerate(order[:150], 1):
            Pres = Pres + terms[n]
            if gap_parts(A - Pres, P - Pres, 5) > 0:
                R = cnt; break
        if R is None:
            print(f"{str([int(band[0]), int(band[1])]):>10} {w2:4.1f} {s:7.2f}   no cert at R<=150")
            continue
        core = order[:R]
        # slot assignment: Bragg comb spacing 2pi/s in u; slot index = round(u/(2pi/s))
        slotidx = set(int(round(np.log(n) / (2 * np.pi / s))) for n in core)
        slots = len(slotidx)
        cps = R / slots
        top = norms[order[0]]
        eff100 = sum(1 for n in order if norms[n] >= top / 100)
        print(f"{str([int(band[0]), int(band[1])]):>10} {w2:4.1f} {s:7.2f} {R:6d} {slots:6d} "
              f"{cps:5.2f} {slots*cps:9.1f} {R*w2*w2:6.1f} {eff100:6d}", flush=True)
print("\nRegistered reading: quadratic derived iff slots ~ 1/w AND cps ~ 1/w; T-growth of R should")
print("localize in one factor. If slots is w-flat or cps is superlinear, the hypothesis fails — report as-is.")

#!/usr/bin/env python3
"""§366: derive C(T) in R ~ C(T)/w^2 — the band dependence of the resonant-core count.

At fixed margin depth (true lmin nearest 1e-8), fixed w, sweep band center T across a ladder;
measure R (norm-greedy exact core, k=5 directional) and form C = R*w^2. Fit C(T) against
candidates: constant | RvM density rho(T)=(1/2pi)log(T/2pi) | log T | sqrt(T).
Registered reading (stated before the fit): if C(T)/rho(T) is FLAT across the ladder, the count
law is R ~ rho(T)/w^2 = (zero density)/w^2 — an a-priori structure law tying the exact core to the
Weyl density (the §335 critical-density picture made quantitative on the STRUCTURE side).
Two widths (w=2, 1.5) to confirm the w^2 scaling holds per-band while C(T) is width-independent.
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

def rho(T):
    return np.log(T / (2 * np.pi)) / (2 * np.pi)

CENTERS = [29, 45, 65, 90, 115, 150, 200, 260, 340, 440]   # band centers T
WID = 30.0
print("C(T) = R*w^2 vs band center T, at margin depth ~1e-8")
print(f"{'w':>4} {'Tc':>5} {'s@1e-8':>7} {'true':>10} {'R':>3} {'coremax':>7} {'C=Rw2':>7} "
      f"{'rho(T)':>7} {'C/rho':>7}")
rows = {2.0: [], 1.5: []}
for w in (2.0, 1.5):
    for Tc in CENTERS:
        b0, b1 = Tc - WID / 2, Tc + WID / 2
        set_grid(b1 + 30.0)
        best = None
        for s in np.arange(6.0, 1.09, -0.05):
            cs = np.arange(b0, b1 + 1e-9, s)
            if len(cs) < 4: continue
            f = L.build_forms(cs, w)
            tm = float(np.linalg.eigvalsh(f['full'])[0])
            if tm < 1e-13: break
            d = abs(np.log10(max(tm, 1e-300)) + 8)
            if best is None or d < best[0]: best = (d, s, cs, f, tm)
        if best is None: continue
        _, s, cs, f, tm = best
        terms = prime_terms_w(cs, w)
        order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
        A, P = f['archpole'], f['prime']
        Pres = np.zeros_like(A); R = None
        for cnt, n in enumerate(order[:200], 1):
            Pres = Pres + terms[n]
            if gap_parts(A - Pres, P - Pres, 5) > 0:
                R = cnt; break
        if R is None: continue
        C = R * w * w
        rT = rho(Tc)
        rows[w].append((Tc, R, C, rT))
        print(f"{w:4.1f} {Tc:5d} {s:7.2f} {tm:+10.2e} {R:3d} {max(order[:R]):7d} {C:7.1f} "
              f"{rT:7.4f} {C/rT:7.1f}", flush=True)

print("\nFits of C(T) (w=2 rows):")
r2 = rows[2.0]
T = np.array([x[0] for x in r2], dtype=float)
C = np.array([x[2] for x in r2], dtype=float)
for name, basis in [("const", np.ones_like(T)),
                    ("rho(T)", rho(T)),
                    ("logT", np.log(T)),
                    ("sqrtT", np.sqrt(T))]:
    coef = float(np.linalg.lstsq(basis[:, None], C[:, None], rcond=None)[0][0][0])
    resid = C - coef * basis
    rel = np.sqrt(np.mean(resid ** 2)) / np.mean(C)
    print(f"  C ~ {coef:8.2f} * {name:7s}   rel-rms residual = {rel:.3f}")
print("\nregistered reading: the basis with smallest rel-rms is the C(T) law; if rho(T) wins and")
print("C/rho is flat across BOTH widths, then R ~ rho(T)/w^2 = (Weyl zero density)/w^2 — the")
print("structure-side count tied a-priori to the zero density. Flat C/rho is the deliverable.")

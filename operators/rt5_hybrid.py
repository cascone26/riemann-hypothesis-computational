#!/usr/bin/env python3
"""§349: THE HYBRID CHAIN — leanest elementary certificate architecture.

  full = A - P = [ (A - Tsmooth) - P_small ] - DT,   P_small = primes n<n0, T = tail n>=n0
  bound = gap_k( (A - Tsmooth) - Pres,  P_small - Pres )  -  eps*||C_unit||_F  -  tailbnd
          with Pres = top-R resonant terms of P_small (norm-greedy), k=5.

Exact-treatment count = R (~40/w^2, polynomial); the rest of the small primes enter only
through the directional 2-block gap bound; the infinite tail enters only through one smooth
integral + one Chebyshev-grade constant. Question: does the hybrid keep FULL REACH (deepest
PSD row) at every w with R polynomial — i.e. does the directional bound survive the smooth
substitution? Falsifier: if the hybrid's reach falls well short of the collective chain's,
the R-term directional structure and the collective tail don't compose, and the lean
architecture fails.
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
XS = np.arange(2, NMAX2 + 1)
DEVREL = np.abs(CUMPSI[2:] - XS) / XS

def tail_bound(w, nprobes):
    uN = np.log(NMAX2)
    I, _ = quad(lambda u: np.exp(u / 2 - u * u * w * w / 4), uN, uN + 200)
    return 2 * 1.04 * I * 4 * w * np.sqrt(np.pi) / (2 * np.pi) * nprobes

def prime_terms_small(cs, w, n0):
    out = {}
    for i, n in enumerate(L.NS):
        if L.NS[i] >= n0: break
        if L.LAMW[i] == 0.0: continue
        u = L.LOGN[i]
        M = np.zeros((len(cs), len(cs)))
        for j in range(len(cs)):
            for k in range(j, len(cs)):
                v = 2.0 * L.LAMW[i] * L.g_pair(cs[j], cs[k], w, u)
                M[j, k] = M[k, j] = v
        out[int(n)] = M
    return out

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

def unit_charge_norm(cs, w, u0):
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
    return float(np.linalg.norm(C))

def gap_bound(M, Q, k=5):
    n = M.shape[0]; k = min(k, n - 1)
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :k], UM[:, k:]
    top = float(np.linalg.eigvalsh(np.diag(lM[:k]) - Ub.T @ Q @ Ub)[0])
    q = float(np.linalg.norm(Ub.T @ Q @ Ur, 2))
    pr = float(np.linalg.eigvalsh(Ur.T @ Q @ Ur)[-1])
    return float(np.linalg.eigvalsh(np.array([[top, -q], [-q, lM[k] - pr]]))[0])

BAND = (14.0, 44.0)
CUT = {2.0: 300, 1.5: 3000, 1.2: 10000}          # winning cuts from §348
print("HYBRID CHAIN (band [14,44]): R exact + directional(k=5) on small primes + smooth tail")
print("+ eps-charge (eps = measured max|psi-x|/x beyond cut; Chebyshev headroom per §348)")
print(f"{'w':>4} {'n0':>6} {'eps':>7} {'s':>5} {'true':>11} {'hybrid':>11} {'R':>3} "
      f"{'coremax':>8} {'status':>10}")
for w in (2.0, 1.5, 1.2):
    n0 = CUT[w]
    u0 = np.log(n0) - 1e-12
    eps = float(np.max(DEVREL[XS >= n0]))
    tb = None
    reach = None; lastrow = None
    for s in np.arange(6.0, 1.19, -0.1):
        cs = np.arange(BAND[0], BAND[1] + 1e-9, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        true_m = float(np.linalg.eigvalsh(f['full'])[0])
        if true_m < 1e-13:
            print(f"{w:4.1f} {n0:6d} {eps:7.4f} {s:5.1f}  (rank limit; hybrid reach s={reach})")
            break
        if tb is None: tb = tail_bound(w, len(cs))
        Tsm = smooth_tail(cs, w, u0)
        charge = eps * unit_charge_norm(cs, w, u0) + tb
        terms = prime_terms_small(cs, w, n0)
        order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
        B = f['archpole'] - Tsm
        Psmall = sum(terms.values())
        Pres = np.zeros_like(B); best = None
        for R, n in enumerate(order[:60], 1):
            Pres = Pres + terms[n]
            b = gap_bound(B - Pres, Psmall - Pres) - charge
            if best is None or b > best[1]: best = (R, b, max(order[:R]))
            if b > 0: break
        R, b, cmx = best
        ok = b > 0
        if ok: reach = s
        print(f"{w:4.1f} {n0:6d} {eps:7.4f} {s:5.1f} {true_m:+11.3e} {b:+11.3e} {R:3d} "
              f"{cmx:8d} {'CERT' if ok else 'fail':>10}", flush=True)
print("\nCompare: collective-full reach (§347) was one grid row above the rank limit at every w,")
print("with 114-1280 exact primes. Hybrid target: same reach with R~40/w^2 exact terms.")

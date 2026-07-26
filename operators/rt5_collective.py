#!/usr/bin/env python3
"""§346: THE COLLECTIVE-TAIL CHAIN — Tao's minor-arc division made operational.

Split P = P_resolved(n < n0, exact) + T(n >= n0, collective). Replace T by its smooth
PNT main term  T_smooth[j,k] = 2 * Int_{u0}^inf e^{u/2} g_pair(c_j,c_k,w,u) du   (d psi -> dx)
and CHARGE the fluctuation  DT = T_exact - T_smooth  by its operator norm:

    lambda_min(A - P)  >=  lambda_min(A - P_res - T_smooth)  -  ||DT||_2  -  tailbnd(>NMAX)

M = A - P_res - T_smooth is fully explicit (arch + pole + finitely many prime terms + one
smooth integral — the §341 minorant machinery applies to all of it); the ONLY arithmetic
beyond the resolved primes is ||DT|| = the psi-fluctuation under the Gaussian-Bragg window —
exactly the object PNT/short-interval theory prices. Today: MEASURE ||DT|| and the reach.

Falsifier: if ||DT|| >> margin at all useful depths for every cut u0, the fluctuation charge
IS the wall (quantified, honest). If reach survives with decent margin, the a-priori program
has a new spine: w-independent resolved core + PNT-priced collective.

Section 2: count-law DECONFOUND — R at FIXED margin depth (true margin nearest 1e-8) vs w
(§345 found edge-R is confounded by where the 0.1 s-grid lands in the plunge).
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

def tail_bound(w, nprobes):
    uN = np.log(NMAX2)
    I, _ = quad(lambda u: np.exp(u / 2 - u * u * w * w / 4), uN, uN + 200)
    return 2 * 1.04 * I * 4 * w * np.sqrt(np.pi) / (2 * np.pi) * nprobes

def prime_block(cs, w, nlo, nhi):
    """exact sum of prime-power terms with nlo <= n < nhi."""
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
    """T_smooth[j,k] = 2 Int_{u0}^inf e^{u/2} g_pair du (PNT main term, d psi -> dx)."""
    n = len(cs)
    M = np.zeros((n, n))
    for j in range(n):
        for k in range(j, n):
            I, _ = quad(lambda u, cj=cs[j], ck=cs[k]: np.exp(u / 2) * L.g_pair(cj, ck, w, u),
                        u0, u0 + 300, limit=400)
            M[j, k] = M[k, j] = 2.0 * I
    return M

BAND = (14.0, 44.0)
print("=" * 100)
print("SECTION 1 — collective-tail chain: cert = lmin(A - Pres - Tsmooth) - ||DT|| - tailbnd")
print("=" * 100)
for w in (2.0, 1.5, 1.2):
    tb = None
    for n0 in (37, 100, 500):
        u0 = np.log(n0) - 1e-12
        reach = None
        print(f"\nw={w}  cut n0={n0}:")
        print(f"{'s':>5} {'true':>11} {'lmin(M)':>11} {'||DT||':>10} {'cert':>11}")
        for s in np.arange(6.0, 1.19, -0.2):
            cs = np.arange(BAND[0], BAND[1] + 1e-9, s)
            if len(cs) < 4: continue
            f = L.build_forms(cs, w)
            true_m = float(np.linalg.eigvalsh(f['full'])[0])
            if true_m < 1e-13:
                print(f"{s:5.1f}  (rank limit)"); break
            if tb is None: tb = tail_bound(w, len(cs))
            Pres = prime_block(cs, w, 2, n0)
            Texact = f['prime'] - Pres
            Tsm = smooth_tail(cs, w, u0)
            DT = Texact - Tsm
            nDT = float(np.linalg.norm(DT, 2))
            M = f['archpole'] - Pres - Tsm
            lmin = float(np.linalg.eigvalsh(M)[0])
            cert = lmin - nDT - tb
            mark = " <-- certified" if cert > 0 else ""
            print(f"{s:5.1f} {true_m:+11.3e} {lmin:+11.3e} {nDT:10.3e} {cert:+11.3e}{mark}",
                  flush=True)
            if cert > 0: reach = s
        print(f"  reach (cert>0): s={reach}   [resolved primes: n<{n0} => "
              f"{int(np.sum((L.NS < n0) & (L.LAMW > 0)))} prime powers]")

print("\n" + "=" * 100)
print("SECTION 2 — count law DECONFOUNDED: R (norm-greedy, k=5) at true margin nearest 1e-8")
print("=" * 100)

def gap_parts(M, Q, k):
    n = M.shape[0]; k = min(k, n - 1)
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :k], UM[:, k:]
    top = float(np.linalg.eigvalsh(np.diag(lM[:k]) - Ub.T @ Q @ Ub)[0])
    q = float(np.linalg.norm(Ub.T @ Q @ Ur, 2))
    pr = float(np.linalg.eigvalsh(Ur.T @ Q @ Ur)[-1])
    return float(np.linalg.eigvalsh(np.array([[top, -q], [-q, lM[k] - pr]]))[0])

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

print(f"{'w':>4} {'s@1e-8':>7} {'true':>11} {'R(k=5)':>7} {'coremax':>8} {'RxW':>6}")
prevR = prevw = None
for w in (3.0, 2.5, 2.0, 1.5, 1.3, 1.2):
    best = None
    for s in np.arange(6.0, 1.19, -0.05):
        cs = np.arange(BAND[0], BAND[1] + 1e-9, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        tm = float(np.linalg.eigvalsh(f['full'])[0])
        if tm < 1e-13: break
        d = abs(np.log10(max(tm, 1e-300)) - (-8))
        if best is None or d < best[0]: best = (d, s, cs, f, tm)
    _, s, cs, f, tm = best
    terms = prime_terms_w(cs, w)
    order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
    A, P = f['archpole'], f['prime']
    Pres = np.zeros_like(A); R = None
    for cnt, n in enumerate(order[:150], 1):
        Pres = Pres + terms[n]
        if gap_parts(A - Pres, P - Pres, 5) > 0:
            R = cnt; break
    cmx = max(order[:R]) if R else -1
    expn = ''
    if prevR and R:
        expn = f"   dlogR/dlog(1/w)={np.log(R / prevR) / np.log(prevw / w):.2f}"
    print(f"{w:4.1f} {s:7.2f} {tm:+11.3e} {str(R):>7} {cmx:8d} {(R or 0) * w:6.1f}{expn}",
          flush=True)
    if R: prevR, prevw = R, w

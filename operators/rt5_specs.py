#!/usr/bin/env python3
"""§343: Roundtable-#5 PRE-REGISTERED SPECS — the terminal deposit, executed.

Registered by the Gate (falsifier form, trend form, exact prime counts — no toy-threshold tests):
  SPEC1  R and marginal buyback d(gap_residual)/dR vs w, band [14,44].
         Kills Tao if buyback decays; confirms/refutes Synth's R=5->8->11->15 (~40/w).
  SPEC2  adaptive-core coupling^2/gap GRID vs w AND band height T.
         Kills uniformity if the ratio climbs with T.
  SPEC3  effective-couplings-per-slot vs w at fixed band (envelope-vs-selectivity race).
         Confirms/refutes Synth's eff10 1->2->2->4->5, product ~1/w^2.
  SPEC4  Dyson-Dirac split: Frobenius Sum|c_k|^2 (honest cost input) vs product
         (redundancy witness); pairwise phase correlations as WITNESS COLUMN ONLY.
         Dyson closes if Frobenius stays ~linear while product grows ~1/w^2.
  SPEC5  per-window Chebyshev psi deviations vs |c_k|, EXACT counts (reach small; VK
         would pass trivially on constants), trend in w only.

All ingredients zero-free (arch+pole+primes; no zeros anywhere in the chain quantities).
Definitions (stated so annotations follow measurements):
  envelope halfwidth in u=log n:  2/w   (weight e^{-u^2 w^2/4} = 1/e at u=2/w)
  eff10/eff100: # prime-power terms with ||P_n||_2 within 10x/100x of the top term
  c_n (off-block coupling): ||Ub^T P_n Ur||_F in the k=3 split of M = A - P_core
  Bragg window for psi: u in [log n - pi/s, log n + pi/s]  (half the slot spacing 2pi/s)
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

# ---------- extend prime cutoff to 3000 (small-w need), rebuild module arrays ----------
NMAX2 = 3000
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
CUMPSI = np.cumsum(lam)          # psi(x) = CUMPSI[floor(x)]

def psi(x):
    xi = int(np.floor(x))
    if xi < 1: return 0.0
    return float(CUMPSI[min(xi, NMAX2)])

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def prime_terms_w(cs, w, thresh=1e-15):
    out = {}
    for i, n in enumerate(L.NS):
        if L.LAMW[i] == 0.0: continue
        u = L.LOGN[i]
        scale = L.LAMW[i] * np.exp(-u * u * w * w / 4)
        if scale < thresh: continue
        M = np.zeros((len(cs), len(cs)))
        for j in range(len(cs)):
            for k in range(j, len(cs)):
                v = 2.0 * L.LAMW[i] * L.g_pair(cs[j], cs[k], w, u)
                M[j, k] = M[k, j] = v
        out[int(n)] = M
    return out

def gap_parts(M, Q, k=3):
    n = M.shape[0]; k = min(k, n - 1)
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :k], UM[:, k:]
    top = float(np.linalg.eigvalsh(np.diag(lM[:k]) - Ub.T @ Q @ Ub)[0])
    q = float(np.linalg.norm(Ub.T @ Q @ Ur, 2))
    pr = float(np.linalg.eigvalsh(Ur.T @ Q @ Ur)[-1])
    bound = float(np.linalg.eigvalsh(np.array([[top, -q], [-q, lM[k] - pr]]))[0])
    return bound, top, q, lM[k] - pr

def find_edge(band, w, Rmax=40, extraR=5):
    """descend s to the last chain-certified spacing; return full state there."""
    last = None
    for s in np.arange(6.0, 1.19, -0.1):
        cs = np.arange(band[0], band[1] + 1e-9, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        true_m = float(np.linalg.eigvalsh(f['full'])[0])
        if true_m < 1e-13:
            return last, ('rank-limit', s)
        terms = prime_terms_w(cs, w)
        order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
        A, P = f['archpole'], f['prime']
        seq = []; certR = None
        Pres = np.zeros_like(A)
        for R in range(1, min(Rmax, len(order)) + 1):
            Pres = Pres + terms[order[R - 1]]
            b, top, q, gden = gap_parts(A - Pres, P - Pres, 3)
            seq.append((R, b, q, gden))
            if b > 0 and certR is None: certR = R
            if certR and R >= certR + extraR: break
        if certR is None:
            return last, ('chain-fail', s)
        last = dict(band=band, w=w, s=s, cs=cs, f=f, terms=terms, order=order,
                    seq=seq, R=certR, true=true_m)
    return last, ('smin', 1.2)

# =====================  SPEC 1 + 3 + 4 + 5 on [14,44]  =====================
BAND1 = (14.0, 44.0)
set_grid(44.0 + 30.0)
WS = [4.0, 3.0, 2.5, 2.0, 1.5]
summary = []
print("=" * 100)
print("SPEC 1/3/4/5 — band [14,44], per w at the chain-certified EDGE spacing")
print("=" * 100)
for w in WS:
    res, stop = find_edge(BAND1, w)
    if res is None:
        print(f"w={w}: no certified edge ({stop})"); continue
    s, R, seq = res['s'], res['R'], res['seq']
    terms, order = res['terms'], res['order']
    # SPEC1: buyback sequence
    print(f"\n--- w={w}  edge s={s:.1f} (stop: {stop[0]} at s={stop[1]:.1f})  "
          f"true margin {res['true']:+.3e}  R_cert={R}  core_max_n={max(order[:R])} ---")
    print("SPEC1 buyback:  R :  b(R)        db(R)=b(R)-b(R-1)   q(off-block)   gap-denom")
    prev = None
    for (Ri, b, q, gden) in seq:
        db = '' if prev is None else f"{b - prev:+12.3e}"
        print(f"   {Ri:3d}: {b:+12.3e} {db:>14s}   {q:12.4e}  {gden:+12.4e}")
        prev = b
    # SPEC3: slots race
    norms = {n: float(np.linalg.norm(terms[n], 2)) for n in order}
    topn = norms[order[0]]
    eff10 = sum(1 for n in order if norms[n] >= topn / 10)
    eff100 = sum(1 for n in order if norms[n] >= topn / 100)
    envhw = 2.0 / w
    print(f"SPEC3 slots: eff10={eff10}  eff100={eff100}  envelope_halfwidth(u)={envhw:.3f}  "
          f"product10={eff10 * envhw:.2f}  product100={eff100 * envhw:.2f}")
    # SPEC4: Frobenius split + phase witness (off-block couplings of the eff100 set)
    Pres = sum(terms[n] for n in order[:R])
    M = res['f']['archpole'] - Pres
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :3], UM[:, 3:]
    sel = [n for n in order if norms[n] >= topn / 100]
    cks, obs = {}, {}
    for n in sel:
        OB = Ub.T @ terms[n] @ Ur
        cks[n] = float(np.linalg.norm(OB))
        obs[n] = OB.ravel()
    frob2 = sum(v * v for v in cks.values())
    cmax, cpair = 0.0, None
    for i in range(len(sel)):
        for j in range(i + 1, len(sel)):
            a, b2 = obs[sel[i]], obs[sel[j]]
            den = np.linalg.norm(a) * np.linalg.norm(b2)
            if den == 0: continue
            c = abs(float(a @ b2)) / den
            if c > cmax: cmax, cpair = c, (sel[i], sel[j])
    print(f"SPEC4 cost split: Sum|c_k|^2 (Frobenius, over eff100 set) = {frob2:.4e}   "
          f"[witness: max pairwise phase corr = {cmax:.3f} at {cpair}]")
    # SPEC5: psi-window deviations for top-5 core terms (exact counts)
    print("SPEC5 psi windows (delta = pi/s in u):  n   ||P_n||     |c_n|      "
          "(psi(b)-psi(a)-(b-a))/(b-a)")
    delta = np.pi / s
    for n in order[:5]:
        u = np.log(n)
        a, b2 = np.exp(u - delta), np.exp(u + delta)
        if b2 <= NMAX2:
            dev = (psi(b2) - psi(a) - (b2 - a)) / (b2 - a)
            devs = f"{dev:+.4f}"
        else:
            devs = "window exceeds sieve"
        print(f"    {n:5d}  {norms[n]:.4e}  {cks.get(n, float('nan')):.4e}   {devs}")
    summary.append((w, s, R, max(order[:R]), eff10, eff100, eff10 * envhw, eff100 * envhw, frob2, cmax))

print("\n" + "=" * 100)
print("SUMMARY vs w (band [14,44], at edge):")
print(f"{'w':>4} {'s_edge':>6} {'R':>3} {'coremax':>7} {'eff10':>5} {'eff100':>6} "
      f"{'prod10':>7} {'prod100':>7} {'Frob^2':>11} {'maxcorr':>7}")
for row in summary:
    print(f"{row[0]:4.1f} {row[1]:6.1f} {row[2]:3d} {row[3]:7d} {row[4]:5d} {row[5]:6d} "
          f"{row[6]:7.2f} {row[7]:7.2f} {row[8]:11.4e} {row[9]:7.3f}")

# =====================  SPEC 2: the w x T grid  =====================
print("\n" + "=" * 100)
print("SPEC 2 — coupling^2/gap GRID (adaptive core, at edge) vs w AND band height")
print("=" * 100)
GBANDS = [(14.0, 44.0), (50.0, 80.0), (100.0, 130.0)]
GWS = [3.0, 2.5, 2.0]
print(f"{'band':>10} {'w':>4} {'s_edge':>6} {'R':>3} {'coremax':>7} {'q':>11} {'gapden':>11} "
      f"{'q^2/gapden':>11} {'bound':>11} {'true':>11}")
for band in GBANDS:
    set_grid(band[1] + 30.0)
    for w in GWS:
        res, stop = find_edge(band, w)
        if res is None:
            print(f"{str(band):>10} {w:4.1f}   no certified edge ({stop})"); continue
        Ri, b, q, gden = res['seq'][res['R'] - 1]
        ratio = q * q / gden if gden > 0 else float('inf')
        print(f"{str([int(band[0]), int(band[1])]):>10} {w:4.1f} {res['s']:6.1f} {res['R']:3d} "
              f"{max(res['order'][:res['R']]):7d} {q:11.4e} {gden:11.4e} {ratio:11.4e} "
              f"{b:+11.3e} {res['true']:+11.3e}", flush=True)
print("\nFalsifier readings: SPEC1 dies if db(R) decays before certification at small w; "
      "SPEC2 kills T-uniformity if q^2/gapden climbs with band; SPEC4 Dyson closes if "
      "Frob^2 ~ linear in 1/w while product ~ 1/w^2; witness corr near 1 would activate "
      "Dirac's theta(gamma)-only escalation gate.")

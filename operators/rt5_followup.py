#!/usr/bin/env python3
"""§345: the three fronts opened by the §344 prints — all zero-free.

A) SELECTION: greedy-by-norm vs greedy-by-BOUND core selection (the §344 non-monotone
   buyback at w=1.5 says norm-ordering is suboptimal). Question: how much R does optimal
   selection save, and does it restore monotone buyback?
B) MERGED SLOTS: cluster u-adjacent prime powers (du < pi/(2s)) into single block units
   (the §344 phase-corr witness hit 0.96 on (4,7),(17,19) — near-parallel off-block
   vectors want joint treatment). Question: does the block chain certify with fewer
   units / monotone buyback?
C) DEEP-W: the count law past w=1.5 (w=1.3,1.2,1.1,1.0), NMAX=20000, explicit
   truncation tail bound printed per row (rows where tail > 0.1x margin are flagged
   TRUNC-LIMITED, not banked). Question: is the w=1.5 exponent acceleration (~3) a
   crossover or a divergence?

Falsifiers: A dies if best-selection R ~= norm R (selection is not the knob);
B dies if block count ~= prime count (merging buys nothing);
C: if d log R / d log(1/w) keeps climbing past ~3, the super-linear law is real
   and 'polynomial wall' on the w-axis is dead; if it flattens, crossover.
"""
import numpy as np
import l2_edge_sweep as L
from scipy.integrate import quad

# ---------- primes to 20000 ----------
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
    """rigorous-shape truncation caveat: entrywise |prime tail beyond NMAX2| bound x nprobes.
    sum_{n>N} Lambda(n) n^-1/2 e^{-u^2 w^2/4} * |gpair-prefactor<=4*w*sqrt(pi)/2pi| * 2,
    psi(x) <= 1.04x."""
    uN = np.log(NMAX2)
    I, _ = quad(lambda u: np.exp(u / 2 - u * u * w * w / 4), uN, uN + 200)
    per_entry = 2 * 1.04 * I * 4 * w * np.sqrt(np.pi) / (2 * np.pi)
    return per_entry * nprobes

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

def gap_parts(M, Q, k=3):
    n = M.shape[0]; k = min(k, n - 1)
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :k], UM[:, k:]
    top = float(np.linalg.eigvalsh(np.diag(lM[:k]) - Ub.T @ Q @ Ub)[0])
    q = float(np.linalg.norm(Ub.T @ Q @ Ur, 2))
    pr = float(np.linalg.eigvalsh(Ur.T @ Q @ Ur)[-1])
    return float(np.linalg.eigvalsh(np.array([[top, -q], [-q, lM[k] - pr]]))[0])

def chain_norm(A, P, units, Rmax=45):
    """units: list of (label, matrix). Greedy by unit norm. Returns (certR, seq, labels, negsteps)."""
    order = sorted(range(len(units)), key=lambda i: -np.linalg.norm(units[i][1], 2))
    seq = []; Pres = np.zeros_like(A); certR = None; neg = 0; prev = None; labs = []
    for cnt, idx in enumerate(order[:Rmax], 1):
        Pres = Pres + units[idx][1]
        labs.append(units[idx][0])
        b = gap_parts(A - Pres, P - Pres)
        if prev is not None and b < prev: neg += 1
        prev = b
        seq.append(b)
        if b > 0:
            certR = cnt; break
    return certR, seq, labs, neg

def chain_best(A, P, units, Rmax=45, cand=45):
    """greedy-by-bound: at each step add the unit that maximizes the new gap bound."""
    order = sorted(range(len(units)), key=lambda i: -np.linalg.norm(units[i][1], 2))[:cand]
    remaining = list(order)
    seq = []; Pres = np.zeros_like(A); certR = None; neg = 0; prev = None; labs = []
    for cnt in range(1, Rmax + 1):
        best = None
        for idx in remaining:
            Pt = Pres + units[idx][1]
            b = gap_parts(A - Pt, P - Pt)
            if best is None or b > best[1]: best = (idx, b)
        idx, b = best
        remaining.remove(idx)
        Pres = Pres + units[idx][1]
        labs.append(units[idx][0])
        if prev is not None and b < prev: neg += 1
        prev = b
        seq.append(b)
        if b > 0:
            certR = cnt; break
    return certR, seq, labs, neg

def cluster_units(terms, s):
    """merge prime powers with du < pi/(2s) into blocks (chain clustering in u)."""
    ns = sorted(terms)
    thresh = np.pi / (2 * s)
    blocks = []; cur = [ns[0]]
    for n in ns[1:]:
        if np.log(n) - np.log(cur[-1]) < thresh: cur.append(n)
        else: blocks.append(cur); cur = [n]
    blocks.append(cur)
    return [('+'.join(map(str, b)), sum(terms[n] for n in b)) for b in blocks]

def edge_s(band, w):
    """last spacing before rank limit (true margin < 1e-13)."""
    prev = None
    for s in np.arange(6.0, 1.19, -0.1):
        cs = np.arange(band[0], band[1] + 1e-9, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        if float(np.linalg.eigvalsh(f['full'])[0]) < 1e-13:
            return prev
        prev = (s, cs, f)
    return prev

BAND = (14.0, 44.0)
print("=" * 100)
print("A + B — selection & merged slots, band [14,44], at each w's edge spacing")
print("=" * 100)
print(f"{'w':>4} {'s':>4} {'true':>10} | {'R_norm':>6} {'neg':>3} {'coremax':>7} | "
      f"{'R_best':>6} {'neg':>3} {'coremax':>7} | {'Rblk_norm':>9} {'neg':>3} | "
      f"{'Rblk_best':>9} {'neg':>3} {'#primes_in_blks':>15}")
resAB = []
for w in (2.5, 2.0, 1.5):
    s, cs, f = edge_s(BAND, w)
    A, P = f['archpole'], f['prime']
    true_m = float(np.linalg.eigvalsh(f['full'])[0])
    terms = prime_terms_w(cs, w)
    prime_units = [(str(n), terms[n]) for n in sorted(terms)]
    Rn, seqn, labn, negn = chain_norm(A, P, prime_units)
    Rb, seqb, labb, negb = chain_best(A, P, prime_units)
    blocks = cluster_units(terms, s)
    Rbn, seqbn, labbn, negbn = chain_norm(A, P, blocks)
    Rbb, seqbb, labbb, negbb = chain_best(A, P, blocks)
    nprim = sum(lab.count('+') + 1 for lab in labbb) if Rbb else -1
    cmx = lambda labs: max(int(x) for lab in labs for x in lab.split('+'))
    print(f"{w:4.1f} {s:4.1f} {true_m:+10.3e} | {str(Rn):>6} {negn:3d} {cmx(labn) if Rn else '-':>7} | "
          f"{str(Rb):>6} {negb:3d} {cmx(labb) if Rb else '-':>7} | {str(Rbn):>9} {negbn:3d} | "
          f"{str(Rbb):>9} {negbb:3d} {nprim:15d}", flush=True)
    resAB.append((w, s, Rn, Rb, Rbn, Rbb, negn, negb, negbn, negbb))
    if w == 1.5:
        print("  w=1.5 best-selection pick order:", ' '.join(labb))
        print("  w=1.5 block pick order (best):  ", ' | '.join(labbb))

print("\n" + "=" * 100)
print("C — DEEP-W: the count law past w=1.5 (best-selection chain, blocks), NMAX=20000")
print("=" * 100)
print(f"{'w':>4} {'s_edge':>6} {'Nz':>3} {'true':>10} {'tailbnd':>9} {'flag':>13} "
      f"{'R_best':>6} {'Rblk_best':>9} {'coremax':>7} {'RxW':>6}")
prevR = None; prevw = None
for w in (1.5, 1.3, 1.2, 1.1, 1.0):
    got = edge_s(BAND, w)
    if got is None:
        print(f"{w:4.1f}  no certified edge"); continue
    s, cs, f = got
    A, P = f['archpole'], f['prime']
    true_m = float(np.linalg.eigvalsh(f['full'])[0])
    tb = tail_bound(w, len(cs))
    flag = "TRUNC-LIMITED" if tb > 0.1 * abs(true_m) else "ok"
    lo, hi = BAND[0] - 2 * w, BAND[1] + 2 * w
    nz = int(np.sum((L.ZEROS >= lo) & (L.ZEROS <= hi)))
    terms = prime_terms_w(cs, w)
    prime_units = [(str(n), terms[n]) for n in sorted(terms)]
    Rb, seqb, labb, negb = chain_best(A, P, prime_units, Rmax=60, cand=60)
    blocks = cluster_units(terms, s)
    Rbb, seqbb, labbb, negbb = chain_best(A, P, blocks, Rmax=60, cand=60)
    cmx = max(int(x) for lab in labb for x in lab.split('+')) if Rb else -1
    expn = ''
    if prevR and Rb:
        expn = f"  dlogR/dlog(1/w)={np.log(Rb / prevR) / np.log(prevw / w):.2f}"
    print(f"{w:4.1f} {s:6.1f} {nz:3d} {true_m:+10.3e} {tb:9.1e} {flag:>13} "
          f"{str(Rb):>6} {str(Rbb):>9} {cmx:7d} {(Rb or 0) * w:6.1f}{expn}", flush=True)
    if Rb: prevR, prevw = Rb, w
print("\nReadings: A) R_best << R_norm => selection is a real knob; B) Rblk < R => merged "
      "slots are the right units; C) exponent flattening => crossover, climbing => the "
      "w-axis wall is super-polynomial and §344's demotion was itself optimistic.")

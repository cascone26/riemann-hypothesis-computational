#!/usr/bin/env python3
"""§345b: WHY does the chain fail at w<=1.3 (R<=60, k=3)? Two suspects, tested head-on:
  (1) the k=3 directional block (banked §339 'k never grows') may finally need to GROW;
  (2) units: bounded-width blocks (du < 0.12, the (17,19)-type correlated pairs) vs per-prime.
Grid: w in {1.5, 1.3, 1.2}, at edge spacing; k in {3,4,5,6,8,10}; norm-greedy core R<=120.
Also prints R needed vs k — if certification returns at k=5 with modest R, the wall statement
becomes 'the DIRECTIONAL block grows as the class completes' (a new, sharper cost law).
"""
import numpy as np
import l2_edge_sweep as L

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

def gap_parts(M, Q, k):
    n = M.shape[0]; k = min(k, n - 1)
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :k], UM[:, k:]
    top = float(np.linalg.eigvalsh(np.diag(lM[:k]) - Ub.T @ Q @ Ub)[0])
    q = float(np.linalg.norm(Ub.T @ Q @ Ur, 2))
    pr = float(np.linalg.eigvalsh(Ur.T @ Q @ Ur)[-1])
    return float(np.linalg.eigvalsh(np.array([[top, -q], [-q, lM[k] - pr]]))[0])

def chain(A, P, units, k, Rmax=120):
    order = sorted(range(len(units)), key=lambda i: -np.linalg.norm(units[i][1], 2))
    Pres = np.zeros_like(A); best = None
    for cnt, idx in enumerate(order[:Rmax], 1):
        Pres = Pres + units[idx][1]
        b = gap_parts(A - Pres, P - Pres, k)
        if best is None or b > best[1]: best = (cnt, b)
        if b > 0:
            return cnt, b
    return None, best[1]

def small_blocks(terms, du=0.12):
    ns = sorted(terms)
    blocks = []; cur = [ns[0]]
    for n in ns[1:]:
        if np.log(n) - np.log(cur[-1]) < du: cur.append(n)
        else: blocks.append(cur); cur = [n]
    blocks.append(cur)
    return [('+'.join(map(str, b)), sum(terms[n] for n in b)) for b in blocks]

def edge_s(band, w):
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
print("k-PROBE at the edge: cert R (or best bound if none) vs directional block size k")
print(f"{'w':>4} {'s':>4} {'n':>3} {'true':>10} {'units':>7} | " +
      " | ".join(f"k={k}" for k in (3, 4, 5, 6, 8, 10)))
for w in (1.5, 1.3, 1.2):
    got = edge_s(BAND, w)
    s, cs, f = got
    A, P = f['archpole'], f['prime']
    true_m = float(np.linalg.eigvalsh(f['full'])[0])
    terms = prime_terms_w(cs, w)
    for units, uname in ((([(str(n), terms[n]) for n in sorted(terms)]), 'prime'),
                         ((small_blocks(terms)), 'blk.12')):
        row = f"{w:4.1f} {s:4.1f} {len(cs):3d} {true_m:+10.3e} {uname:>7} | "
        for k in (3, 4, 5, 6, 8, 10):
            R, b = chain(A, P, units, k)
            row += (f"R={R:<4d}" if R else f"x{b:+.0e}") + " | "
        print(row, flush=True)
    # sanity: max k possible is n-1; note largest block content for blk.12
    bl = small_blocks(terms)
    big = max(bl, key=lambda t: t[0].count('+'))
    print(f"      blk.12: {len(bl)} units; largest block has {big[0].count('+')+1} members "
          f"(starts {big[0].split('+')[0]})", flush=True)

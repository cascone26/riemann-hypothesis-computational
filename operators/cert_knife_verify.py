#!/usr/bin/env python3
"""§832-VERIFY: is the 4W=2.2 margin collapse (0.103→5.6e-5) REAL prime-side thinning or a bad-span
conditioning artifact? Re-certify at 4W∈{1.60,1.90,2.20} sweeping N, reporting the certified Weil
pencil-min AND G-conditioning (λmin(G), cond(G), float midpoint proxy). Decisive test: at 4W=2.2 does
the margin stay collapsed even at a conservative, WELL-conditioned N → real thinning; or recover when
G is healthy → conditioning. Reuses §830/§832 arb entry functions."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec = 250

def certmargin(fourW, CEN, rc=400):
    ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
    cs = [arb(c) for c in CEN]; n = len(cs)
    Q = [[None]*n for _ in range(n)]; G = [[None]*n for _ in range(n)]; maxrad = 0.0
    Gf = np.zeros((n, n))
    for j in range(n):
        for k in range(j, n):
            Q[j][k] = Q[k][j] = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(rc)) - ac.prime_entry(cs[j], cs[k])
            G[j][k] = G[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))
            Gf[j, k] = Gf[k, j] = float(G[j][k].mid())
            maxrad = max(maxrad, float(Q[j][k].rad()))
    gw = np.linalg.eigvalsh(Gf); lminG = float(gw[0]); condG = float(gw[-1]/max(gw[0], 1e-300))
    def ok(lam):
        d = [None]*n; L = [[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd = Q[i][i] - lam*G[i][i]
            for kk in range(i): sd -= L[i][kk]*L[i][kk]*d[kk]
            if not (sd > 0): return False
            d[i] = sd
            for j2 in range(i+1, n):
                v = Q[j2][i] - lam*G[j2][i]
                for kk in range(i): v -= L[j2][kk]*L[i][kk]*d[kk]
                L[j2][i] = v/d[i]
        return True
    if not ok(arb(0)): return None, maxrad, lminG, condG
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    return lo, maxrad, lminG, condG

LOGN = np.log(np.arange(2, 300))
print("=" * 104)
print("§832-VERIFY — is the 4W=2.2 margin collapse real prime-side thinning or bad-span conditioning?")
print("=" * 104)
print(f"{'4W':>6} {'N':>4} {'N/Ncap':>7} {'nprime':>7} {'cert margin':>14} {'maxrad':>9} {'λmin(G)':>10} {'cond(G)':>10}  verdict")
for fourW in [1.60, 1.90, 2.20]:
    Ncap = 15*fourW/np.pi; npr = int((LOGN < fourW).sum())
    # sweep N from conservative (well below cap) up to ~0.8 cap
    Ns = sorted(set(max(4, round(f*Ncap)) for f in [0.55, 0.65, 0.76]))
    print("-"*104)
    for N in Ns:
        CEN = list(np.linspace(14, 44, N))
        m, mr, lg, cg = certmargin(fourW, CEN)
        val = f"{m:.3e}" if m is not None else "swallowed"
        health = "G-healthy" if (lg > 0.5 and cg < 50) else ("G-loading" if cg < 1e4 else "G-DEGEN")
        print(f"{fourW:>6.2f} {N:>4} {N/Ncap:>7.2f} {npr:>7} {val:>14} {mr:>9.1e} {lg:>10.3f} {cg:>10.1f}  {health}", flush=True)
print("=" * 104)
print("READ: if 4W=2.2 margin stays small (~1e-4) even at conservative N with G-healthy → REAL prime thinning.")
print("      if margin recovers to ~0.1 at well-conditioned N → the N=8 collapse was conditioning.")
print("=" * 104)

#!/usr/bin/env python3
"""PASS C (Euler) — c vs window length L_N in Nyquist units, at fixed 4W=1.32.
Pass A/B (frame_const_controlled.out) showed c is independent of W and of density excess at
fixed L_N=4, fluctuating only with zero-grid phase. So the remaining variable is L_N itself.
n probes at exact Nyquist spacing from floor t*+2; n = 4,6,8,10,12,16.
Registered prediction: c decays with L_N toward a positive floor governed by the excess
(Beurling: excess>0 suffices for an infinite-window frame); the decay profile IS the
finite-defect content the frame lemma / horizon computation needs.
"""
import sys, os, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac

ZC = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros500.npy'))

def cert(fourW, centers, gh=800):
    ctx.prec = 220
    ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
    Wa = ac.W; cs = [arb(c) for c in centers]; n = len(cs)
    def fej_arb(c, r):
        out = arb(0)
        for cc in (c, -c):
            y = Wa*(r-cc)
            out = out + (arb(1).union((y.sin()/y)**2) if y.contains(arb(0)) else (y.sin()/y)**2)
        return out
    S = [[arb(0)]*n for _ in range(n)]
    for g in ZC:
        if g > gh: break
        gb = arb(float(g), 1e-9)
        F = [fej_arb(cs[j], gb) for j in range(n)]
        for j in range(n):
            for k in range(j, n): S[j][k] = S[j][k] + F[j]*F[k]
    for j in range(n):
        for k in range(j, n): S[k][j] = S[j][k]
    GM = [[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j, n): GM[j][k] = GM[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))
    def ok(lam):
        d = [None]*n; L = [[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd = S[i][i] - lam*GM[i][i]
            for k in range(i): sd -= L[i][k]*L[i][k]*d[k]
            if not (sd > 0): return False
            d[i] = sd
            for j2 in range(i+1, n):
                v = S[j2][i] - lam*GM[j2][i]
                for k in range(i): v -= L[j2][k]*L[i][k]*d[k]
                L[j2][i] = v/d[i]
        return True
    if not ok(arb(0)): return None
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    return lo

fourW = 1.32
ts = 2*math.pi*math.exp(fourW)
sp = 2*math.pi/fourW
floor = ts + 2.0
print(f"4W={fourW}, floor={floor:.2f}, Nyquist sp={sp:.2f}; c vs n probes (L_N = n):")
for n in (4, 6, 8, 10, 12, 16):
    centers = [floor + i*sp for i in range(n)]
    c = cert(fourW, centers)
    top = centers[-1]+sp/2
    print(f"  n={n:>2}  window [{floor-sp/2:.1f},{top:.1f}]  "
          + (f"c = {c:.6e}" if c is not None else "S NOT PD"), flush=True)
print("EXIT:0")

#!/usr/bin/env python3
"""Extend §796: certified frame constants at more (4W, window) anchors. Capacity-matched probes
(spacing ~ 2pi/4W), windows above t*(W). Banks as §796-addendum."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac

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
    ZC = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros500.npy'))
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
    if not ok(arb(0)):
        print(f"4W={fourW} centers {centers}: S NOT PD (over capacity?)", flush=True); return
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    print(f"4W={fourW} centers {centers}: CERTIFIED c = {lo:.6e}", flush=True)

import math
for fourW, cmax in ((1.00, 44), (1.32, 44), (1.61, 60), (2.00, 72)):
    ts = 2*math.pi*math.exp(fourW)
    sp = max(4.0, 2*math.pi/fourW*1.0)
    centers = []
    c = math.ceil(ts+2)
    while c <= cmax: centers.append(round(c,1)); c += sp
    if len(centers) >= 3: cert(fourW, centers)
    else: print(f"4W={fourW}: window above t*={ts:.1f} too small for 3 probes at cmax={cmax}")

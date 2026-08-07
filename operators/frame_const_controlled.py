#!/usr/bin/env python3
"""CONTROLLED c(W) scan (Euler, round 3) — §796-curve confound removed.
frame_const_curve.py confounds two variables: sp = max(4.0, 2pi/4W) puts the 2.00 anchor at
fill 0.785 of Nyquist (sub-capacity => inflated c) while 1.00/1.32 sit at fill 1.000; and
ceil(t*+2) rounds the floor by a varying phase fraction. This scan controls both:
  n = 4 probes ALWAYS, sp = 2pi/4W EXACT (fill 1.000), floor = t*(W) + 2.0 EXACT.
Pass B: floor shifted by +sp/2 at three widths — isolates grid-phase fluctuation at fixed W.
Registered predictions (before run): H1 (c ~ mean density excess) => c falls ~5-6x from
4W=0.95 to 2.25. H2 (c = O(1), phase-fluctuating) => flat with O(1) scatter, pass B moves c
as much as moving W does. Same certified pencil machinery as §796 (bisected shifted LDL^T,
S - c*G, hypothesis: ordinate balls only).
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
    if not ok(arb(0)):
        return None
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    return lo

def diag(fourW, centers):
    sp = centers[1]-centers[0]
    lo, hi = centers[0]-sp/2, centers[-1]+sp/2
    zin = [g for g in ZC if lo <= g <= hi]
    exc = np.mean([math.log(t/(2*math.pi))-fourW for t in np.linspace(lo, hi, 200)])
    cap = (hi-lo)*fourW/(2*math.pi)
    return len(zin), cap, exc

print("PASS A — controlled: n=4, sp=2pi/4W exact, floor=t*+2 exact")
print(f"{'4W':>5} {'floor':>7} {'sp':>6} {'nz':>3} {'cap':>5} {'excess':>7} {'c':>12}")
for fourW in (0.95, 1.00, 1.10, 1.20, 1.32, 1.45, 1.61, 1.80, 2.00, 2.25):
    ts = 2*math.pi*math.exp(fourW)
    sp = 2*math.pi/fourW
    floor = ts + 2.0
    centers = [floor + i*sp for i in range(4)]
    c = cert(fourW, centers)
    nz, cap, exc = diag(fourW, centers)
    print(f"{fourW:>5} {floor:>7.2f} {sp:>6.2f} {nz:>3} {cap:>5.2f} {exc:>7.3f} "
          + (f"{c:>12.6e}" if c is not None else "   S NOT PD"), flush=True)

print("\nPASS B — phase test: same widths, floor shifted +sp/2")
for fourW in (1.00, 1.32, 2.00):
    ts = 2*math.pi*math.exp(fourW)
    sp = 2*math.pi/fourW
    floor = ts + 2.0 + sp/2
    centers = [floor + i*sp for i in range(4)]
    c = cert(fourW, centers)
    nz, cap, exc = diag(fourW, centers)
    print(f"{fourW:>5} {floor:>7.2f} {sp:>6.2f} {nz:>3} {cap:>5.2f} {exc:>7.3f} "
          + (f"{c:>12.6e}" if c is not None else "   S NOT PD"), flush=True)
print("EXIT:0")

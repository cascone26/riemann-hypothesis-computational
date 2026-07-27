#!/usr/bin/env python3
"""§379b: stiffness at FIXED margin depth (kills the §379 depth-aliasing).
For each band, bisect s so true lmin(full) ~ target (1e-6). Then bend Lambda(2)'s powers up until
lmin crosses 0; report b*-1 = the Lambda(2) cushion at matched depth. Band-flat b* => uniform
arithmetic stiffness at fixed depth; decaying => cushion shrinks with T even at matched margin.
Also report the cushion normalized by the n=2 term-norm (removes the trivial magnitude scaling)."""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def scale2(b):
    s = np.ones_like(L.LAMW)
    for i, n in enumerate(L.NS):
        m = n
        while m % 2 == 0: m //= 2
        if m == 1: s[i] = b
    return s

w = 3.0
TARGET = 1e-6
BANDS = [(14.0, 44.0), (30.0, 60.0), (50.0, 80.0), (75.0, 105.0), (100.0, 130.0), (135.0, 165.0)]
print(f"§379b stiffness at FIXED depth lmin0~{TARGET:.0e} (w=3): Lambda(2) cushion b*-1 vs band")
print(f"{'band':>10} {'s':>5} {'lmin0':>10} {'b*-1':>8} {'n2termnorm':>11} {'cushion*norm':>12}")
for band in BANDS:
    b0, b1 = band
    set_grid(b1 + 30.0)
    # bisect s in [1.2, 6.5] to hit lmin ~ TARGET (lmin decreasing in #probes => increasing as s grows)
    def lmin_s(s):
        cs = np.arange(b0, b1 + 1e-9, s)
        if len(cs) < 4: return None
        return float(np.linalg.eigvalsh(L.build_forms(cs, w)['full'])[0]), cs
    lo, hi = 1.5, 6.0
    ss = None
    for _ in range(45):
        mid = (lo + hi) / 2
        r = lmin_s(mid)
        if r is None or r[0] is None:
            lo = mid; continue
        lm = r[0]
        if lm > TARGET: hi = mid      # larger s -> larger margin; too big -> lower s
        else: lo = mid
        ss = mid
    s = (lo + hi) / 2
    lm, cs = lmin_s(s)
    if lm is None or lm <= 0:
        print(f"{str([int(b0),int(b1)]):>10} depth miss"); continue
    # bend Lambda(2)
    def lmin_b(b):
        return float(np.linalg.eigvalsh(L.build_forms(cs, w, lam_scale=scale2(b))['full'])[0])
    blo, bhi = 1.0, 1.0
    ok = True
    for _ in range(80):
        bhi *= 1.03
        if lmin_b(bhi) <= 0: break
    else: ok = False
    if ok:
        for _ in range(50):
            bm = (blo + bhi) / 2
            if lmin_b(bm) > 0: blo = bm
            else: bhi = bm
        bstar = (blo + bhi) / 2 - 1.0
    else:
        bstar = None
    # n=2 term norm on this config
    u2 = np.log(2); n2 = np.zeros((len(cs), len(cs)))
    for j in range(len(cs)):
        for k in range(j, len(cs)):
            n2[j, k] = n2[k, j] = 2.0 * (np.log(2) / np.sqrt(2)) * L.g_pair(cs[j], cs[k], w, u2)
    tn = float(np.linalg.norm(n2, 2))
    fmt = lambda x: f"{x:+.4f}" if x is not None else " none "
    cn = f"{bstar*tn:.4e}" if bstar is not None else "  -  "
    print(f"{str([int(b0),int(b1)]):>10} {s:5.2f} {lm:10.2e} {fmt(bstar)} {tn:11.4e} {cn:>12}", flush=True)
print("\nreading: b*-1 band-flat at matched depth => uniform Lambda(2) stiffness (RH-robustness); the")
print("magnitude-normalized cushion (b*-1)*termnorm removes the trivial scaling and is the intrinsic cushion.")

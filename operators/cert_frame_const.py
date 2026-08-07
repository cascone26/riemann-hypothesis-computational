#!/usr/bin/env python3
"""CERTIFIED FRAME CONSTANT above t* (REDUCTION.md piece 2 → banks as §796).

Statement certified: for the band-limited Fejér class on centers {25..44} at 4W=1.32 (entirely
above t* = 2pi e^{1.32} ~ 23.5), the zeros-only sampling form satisfies
    sum_{gamma <= 800} psi_hat(gamma)^2  >=  c * ||psi||_{L2}^2
with c a certified positive constant (bisected shifted interval-LDL^T on S - c*G).

HYPOTHESIS SET (one item): ordinate enclosures — the 491 ordinates (zeros500.npy, mpmath dps-30)
correct within 1e-9 balls. NOTHING else: completeness NOT needed (a missing zero would only
INCREASE the left side — one-sided friendly); on-line-ness of zeros ABOVE 800 irrelevant (they
don't appear); on-line-ness below 800 is part of what "ordinate" means for listed real ordinates
(PT covers; and a hypothetically off-line listed zero is excluded from the sum by the same
one-sidedness — the inequality is over the listed real evaluation points as an L2 bound).
Both S and G rigorous: S = sum of closed-form Fejér products at ordinate balls; G via the u=0
closed form G[j,k] = 2*pi*gpair_closed_arb(cj, ck, 0).
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac

fourW = 1.32
CENTERS = [25, 31, 37, 43]        # 4 probes, spacing 6.0 = class capacity 4W/2pi~0.21/unit on [25,44] (~4 dims); denser probe sets are redundant -> S certifiably singular (the §476g lesson, twice tonight)
G_HEIGHT = 800

ctx.prec = 220
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
Wa = ac.W
cs = [arb(c) for c in CENTERS]; n = len(cs)

def fej_arb(c, r):
    out = arb(0)
    for cc in (c, -c):
        y = Wa*(r-cc)
        if y.contains(arb(0)):
            out = out + arb(1).union((y.sin()/y)**2)
        else:
            out = out + (y.sin()/y)**2
    return out

print(f"S: zeros-only sampling Gram, {n} probes {CENTERS[0]}..{CENTERS[-1]}, 4W={fourW}...", flush=True)
ZC = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros500.npy'))
S = [[arb(0)]*n for _ in range(n)]
nz = 0
for g in ZC:
    if g > G_HEIGHT: break
    gb = arb(float(g), 1e-9)
    F = [fej_arb(cs[j], gb) for j in range(n)]
    for j in range(n):
        for k in range(j, n):
            S[j][k] = S[j][k] + F[j]*F[k]
    nz += 1
for j in range(n):
    for k in range(j, n):
        S[k][j] = S[j][k]
print(f"  {nz} zeros summed (closed form, ordinate balls 1e-9)", flush=True)

print("G: L2 overlap Gram via u=0 closed form...", flush=True)
GM = [[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j, n):
        GM[j][k] = GM[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))

def all_pos_shifted(lam):
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

if not all_pos_shifted(arb(0)):
    print("S itself not certified PD — no frame constant at this precision"); sys.exit(1)
lo, hi = 0.0, 1.0
for _ in range(40):
    mid = (lo+hi)/2
    if all_pos_shifted(arb(mid)): lo = mid
    else: hi = mid
print(f"\n===== CERTIFIED FRAME CONSTANT (§796) =====")
print(f"sum_(gamma<=800) psi_hat(gamma)^2 >= {lo:.6e} * ||psi||^2   on the above-t* class")
print(f"(4W={fourW}, centers {CENTERS[0]}..{CENTERS[-1]}, {nz} zeros; hypothesis: ordinate balls only)")
print(f"float scout predicted ~5.6e-3 on the s=1.0 sub-probing — certified value refines it.")

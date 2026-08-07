#!/usr/bin/env python3
"""RIEMANN8 item 1a — Landau/jump-energy ratios of the §476c certified witness, certified balls.
Identity (§473/§476d): E_n(psi) = w_n * int (psi(u)-psi(u+log n))^2 du = 2 w_n (||psi||^2 - ghat(log n)),
so  E_n/(2 w_n ||psi||^2) = 1 - r_n,  r_n = ghat_vv(log n)/ghat_vv(0).
Required for full(v) >= 0 given full_no3(v) < 0:  E_3-ratio >= 1 + |Q_no3(v)|/(2 w_3 ||v||^2).
All three ghat values are CLOSED FORM (arb_fejer_cert2.gpair_closed_arb) — no integrals, no zeros.
Witness: identical deterministic reconstruction to arb_rung_certA2.py (float min-eig of no-3 form).
Q_no3(v): banked certified enclosure §476c: [-8.45e-5, -0.8e-5].
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import l2_fejer as fj

fourW = 1.32; S = 2.0
ctx.prec = 200

# --- witness, exactly as certA2 ---
fj.W = fourW/4.0
cs_f = fj.probes((14.0, 44.0), S)
scale = np.ones_like(fj.LAMW); scale[np.isin(fj.NS, [3, 9, 27, 81])] = 0.0
fmat = fj.build(cs_f, lam_scale=scale)
FU = 0.5*(fmat['full']+fmat['full'].T)
ev, V = np.linalg.eigh(FU)
print(f"float: n={len(cs_f)} probes (s={S}), min eig(no-3) = {ev[0]:.4e}  (§476c: -8.52e-5)")
v = V[:, 0]

# --- arb setup, as certA2 ---
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
cs = [arb(float(c)) for c in cs_f]
va = [arb(float(x)) for x in v]
n = len(cs)

def ghat(u):
    g = arb(0)
    for j in range(n):
        for k in range(n):
            g = g + va[j]*va[k]*ac.gpair_closed_arb(cs[j], cs[k], u)
    return g

g0 = ghat(arb(0))
g2 = ghat(arb(2).log())
g3 = ghat(arb(3).log())
print(f"ghat_vv(0)    = {g0}")
print(f"ghat_vv(log2) = {g2}")
print(f"ghat_vv(log3) = {g3}")

r2 = g2/g0; r3 = g3/g0
R2 = 1 - r2; R3 = 1 - r3
print(f"\nr_log2 = {r2}   (banked float §476f: +0.413)")
print(f"r_log3 = {r3}   (banked float §476e.3/f: -0.0467)")
print(f"\nE2-ratio = 1 - r_log2 = {R2}")
print(f"E3-ratio = 1 - r_log3 = {R3}")

# required threshold from banked certified Q_no3 enclosure (§476c)
w3 = arb(3).log()/arb(3).sqrt()
Qmid = (-8.45e-5 + -0.8e-5)/2; Qrad = (8.45e-5 - 0.8e-5)/2
Qno3 = arb(Qmid, Qrad)
req = 1 + (-Qno3)/(2*w3*g0)
print(f"\nrequired E3-ratio (= 1 + |Q_no3|/(2 w3 ||v||^2), Q_no3 banked §476c) = {req}")
slack = R3 - req
print(f"slack (achieved - required) = {slack}")
print(f"\nfull(v)/(2 w3 ||v||^2) consistency check (= slack): "
      f"{(Qno3 + 2*w3*(g0 - g3))/(2*w3*g0)}")

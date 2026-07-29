#!/usr/bin/env python3
"""§394 (angle 7): Euler-Lagrange criticality — is true Lambda margin-critical, and what is the
gradient dlmin/dLambda(n) (the E-L residual / dual object)? At a near-critical config, perturb each
prime's Lambda and measure dlmin. If true Lambda is extremal, the gradient pattern IS the Lagrange
multiplier of the pair-constraints = candidate new dual arithmetic object."""
import numpy as np
import l2_edge_sweep as L
w = 3.0
def lmin_scaled(cs, scale):
    return float(np.linalg.eigvalsh(L.build_forms(cs, w, lam_scale=scale)['full'])[0])
# near-critical row on [14,44]
for s in np.arange(6.0, 2.9, -0.1):
    cs = np.arange(14.0, 44.0+1e-9, s)
    l0 = float(np.linalg.eigvalsh(L.build_forms(cs, w)['full'])[0])
    if l0 < 1e-4: break
print(f"config s={s:.1f} n={len(cs)} lmin0={l0:.3e}. Gradient dlmin/dlog Lambda(n) per prime power:")
eps = 1e-4
grads = {}
for i, n in enumerate(L.NS):
    if L.LAMW[i] == 0 or L.LAMW[i]*np.exp(-(np.log(n)**2)*w*w/4) < 1e-13: continue
    sc = np.ones_like(L.LAMW); sc[i] = 1+eps
    lp = lmin_scaled(cs, sc); sc[i] = 1-eps; lm = lmin_scaled(cs, sc)
    grads[int(n)] = (lp - lm)/(2*eps)
tot = sum(grads.values())
print(f"total dlmin/dscale (all primes together) approx = {tot:+.4e}")
for n in sorted(grads, key=lambda k:-abs(grads[k]))[:12]:
    print(f"  n={n:4d}: dlmin/dlogLam = {grads[n]:+.4e}  (Lam(n)/sqrt(n) weight {L.LAMW[list(L.NS).index(n)]:.3e})")
print("reading: gradient concentrated + structured = extremal/dual object; if it just tracks term")
print("magnitude (n=2 biggest), it's trivial scaling. Check: is sign uniform (all lower margin) = ")
print("true Lambda sits where ANY increase hurts = margin-critical (two-floor prediction).")

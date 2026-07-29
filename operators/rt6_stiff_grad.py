#!/usr/bin/env python3
"""§412b: the INTRINSIC single-prime stiffness = gradient dλmin/dε_p at ε=1, at a fixed well-conditioned s.
Separates two confounds the threshold (§412) mixed:
  (A) sensitivity per prime = |dλmin/dε_p|  (linear response, threshold-distance-independent)
  (B) #powers of p in the prime-power range NS (low primes have more powers -> more perturbed terms)
Report both, and the PER-POWER-NORMALIZED sensitivity |dλmin/dε_p| / (#powers of p). If the normalized
sensitivity is ~flat across p, the p-grading is just 'more low-prime powers'; if it still falls with p,
the low primes are intrinsically the resonant/mean direction (matches §411)."""
import numpy as np
import l2_edge_sweep as L

BAND = (14.0, 44.0); W = 3.0
PRIMES = (2, 3, 5, 7, 11, 13)

def scale_for(p, eps):
    s = np.ones_like(L.LAMW)
    for i, n in enumerate(L.NS):
        m = n
        while m % p == 0: m //= p
        if m == 1: s[i] = eps
    return s

def npow(p):  # #pure powers of p present in NS
    c = 0
    for n in L.NS:
        m = n
        while m % p == 0: m //= p
        if m == 1: c += 1
    return c

def lmin(cs, w, scale):
    return float(np.linalg.eigvalsh(L.build_forms(cs, w, lam_scale=scale)['full'])[0])

print(f"§412b intrinsic stiffness gradient. band {BAND} w={W}. (2C2-1={2*0.6601618158-1:.3f})")
print(f"{'s':>5} {'n':>3} {'rcond':>8} {'l0':>10}  " + "  ".join(f"p={p}" for p in PRIMES))
h = 1e-4
for s in (6.0, 5.5, 5.0, 4.5):
    cs = L.probes(BAND, s)
    if len(cs) < 4: continue
    f0 = L.build_forms(cs, W); ev = np.linalg.eigvalsh(f0['overlap'])
    rc = ev[-1] / ev[0]; l0 = float(np.linalg.eigvalsh(f0['full'])[0])
    grads = []
    for p in PRIMES:
        lp = lmin(cs, W, scale_for(p, 1 + h)); lm = lmin(cs, W, scale_for(p, 1 - h))
        grads.append((lp - lm) / (2 * h))
    print(f"{s:5.1f} {len(cs):3d} {rc:8.1e} {l0:+10.2e}  " +
          "  ".join(f"{g:+.3f}" for g in grads), flush=True)

print("\n#pure powers of p in NS:", {p: npow(p) for p in PRIMES})
print("per-power normalized |grad| (at s=6.0):")
cs = L.probes(BAND, 6.0)
for p in PRIMES:
    lp = lmin(cs, W, scale_for(p, 1 + h)); lm = lmin(cs, W, scale_for(p, 1 - h))
    g = (lp - lm) / (2 * h)
    print(f"  p={p:2d}: grad={g:+.4f}  npow={npow(p)}  |grad|/npow={abs(g)/npow(p):.4f}  "
          f"Lam(p)/sqrt(p)={np.log(p)/np.sqrt(p):.4f}")
print("""
reading: dλmin/dε_p is the response of the Weil margin to inflating prime p's whole power-tower.
  NEGATIVE grad = inflating p LOWERS the margin (p is a positivity-eroding direction).
  If |grad| falls with p even AFTER dividing by #powers, the low primes are intrinsically the soft
  (resonant/mean) direction -> the sign half is governed by the SMALL primes, not a flat arithmetic floor.
""")

#!/usr/bin/env python3
"""RUNG BUILD 1 (§475 primary lane): the L>=log3 rung, empirical margin map + razor-tuning stiffness.
Fixes the §475d caveat: margin = min over UNIT-L2-NORM test functions in the probe span, computed as the
generalized min-eig of (full, G) on G's well-conditioned subspace (projecting out the probe family's own
near-null directions, which faked machine-zero margins in support_ladder_sweep.py).

Rungs (products F_jF_k bandlimited to |u|<4W => primes log n < 4W):
  4W < log2          : archimedean only  -- PROVEN positive (CC 2006.13771)   [GATE 1]
  log2 <= 4W < log3  : + p=2 only        -- PROVEN positive (CC 2106.01715)   [GATE 2]
  4W >= log3         : + p=3             -- OPEN (the rung)
Expected: margins >= 0 everywhere (zeros are on the line as far as sampled); the INFORMATION is the margin
CURVE (where it collapses toward the a=1 razor edge) and the STIFFNESS: CC report positivity at L=log3
survives only if the p=2 weight is within ~1e-3 of its true value. We test that in our independent
instrument: scale the p=2 Euler weight by (1+eps) (all powers of 2) and find the eps-interval where the
margin stays >= 0, as a function of 4W. Razor-tuning = that interval crashing as 4W crosses log3.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import l2_fejer as fj

BAND = (14.0, 44.0); SP = 1.5
GTOL = 1e-8      # relative Gram eigenvalue cutoff for the well-conditioned subspace

def margin(fourW, eps2=0.0):
    """Generalized min-eig of (full, G) on G's well-conditioned subspace.
    eps2: relative perturbation of the prime-2 Euler weight (n=2,4,8,... scaled by 1+eps2)."""
    fj.W = fourW/4.0
    cs = fj.probes(BAND, SP)
    scale = np.ones_like(fj.LAMW)
    if eps2 != 0.0:
        pow2 = np.isin(fj.NS, [2,4,8,16,32,64,128])
        scale[pow2] = 1.0 + eps2
    f = fj.build(cs, lam_scale=scale)
    FU = 0.5*(f['full']+f['full'].T); G = 0.5*(f['overlap']+f['overlap'].T)
    ge, gv = np.linalg.eigh(G)
    keep = ge > GTOL*ge[-1]
    B = gv[:, keep] / np.sqrt(ge[keep])[None, :]     # G-orthonormal basis of the good subspace
    M = B.T @ FU @ B
    return float(np.linalg.eigvalsh(0.5*(M+M.T))[0]), int(keep.sum()), len(cs)

log2, log3, log4, log5 = np.log(2), np.log(3), np.log(4), np.log(5)
print(f"log2={log2:.4f} log3={log3:.4f} log4={log4:.4f} log5={log5:.4f}\n")
print("PART 1 — margin curve (unit-L2-normalized Weil-form minimum) vs 4W:")
print(f"{'4W':>6} {'rung':>16} {'margin':>13} {'dim':>4}")
for fourW, rung in [(0.55,'arch GATE1'),(0.66,'arch GATE1'),(0.80,'p2 GATE2'),(0.95,'p2 GATE2'),
                    (1.07,'p2 GATE2'),(1.12,'+p3 OPEN'),(1.20,'+p3 OPEN'),(1.32,'+p3 OPEN'),
                    (1.42,'+p4 OPEN'),(1.65,'+p5 OPEN')]:
    m, d, n = margin(fourW)
    print(f"{fourW:>6.2f} {rung:>16} {m:>13.4e} {d:>4}", flush=True)

print("\nPART 2 — razor-tuning: eps-interval of p=2-weight perturbation keeping margin>=0, vs 4W:")
print("(CC report ~1e-3 tolerance at L=log3; our independent check)")
print(f"{'4W':>6} {'m(eps=0)':>12} {'m(+1e-2)':>12} {'m(-1e-2)':>12} {'m(+1e-3)':>12} {'m(-1e-3)':>12}")
for fourW in [0.95, 1.07, 1.12, 1.20, 1.32]:
    row = [margin(fourW, e)[0] for e in [0.0, 1e-2, -1e-2, 1e-3, -1e-3]]
    print(f"{fourW:>6.2f} " + " ".join(f"{x:>12.4e}" for x in row), flush=True)
print("\nReading: GATES 1-2 must show margin>=0 (proven rungs). The open-rung information is (i) the margin")
print("collapse profile crossing log3=1.0986 and (ii) whether the tolerated eps-window crashes there")
print("(= CC's razor-tuning, = our §457/§470, reproduced in an independent instrument).")

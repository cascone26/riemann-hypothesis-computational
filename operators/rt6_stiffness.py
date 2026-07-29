#!/usr/bin/env python3
"""§412 (ATTACK, positivity-bound family): QUANTIFY the critical single-prime stiffness.
§411 showed: TRUE Lambda holds B_analytic PSD to the rank edge; a SYSTEMATIC inflation of Lambda(2)
breaks it EARLY (before the edge); random jitter does not. Now MEASURE the threshold:
  - For each prime p (2,3,5,7), sweep the multiplicative inflation factor eps on Lambda(p) only,
    eps in {1.00 ... 1.20}, and find eps* = the smallest inflation at which lmin(B_analytic) first
    goes negative at a FIXED sub-edge density s (well-conditioned, rcond not tiny).
  - Report eps*-1 (the fractional stiffness margin) per prime, and the WEIGHTED margin
    (eps*-1)*Lambda(p)/sqrt(p) (the actual coefficient perturbation size).
  - Compare against the §385 two-floor relative bound 2*C2-1 ~ 0.3238 (twin-prime constant floor).
QUESTION: is the true form's stiffness (a) prime-uniform, (b) matched to 2C2-1, or (c) p-dependent
(strongest at p=2, the resonant term)? A prime-uniform floor matching 2C2-1 would tie the analytic
positivity margin directly to the arithmetic Parity-Bar constant = the mechanism made quantitative.
"""
import numpy as np
import l2_edge_sweep as L

C2 = 0.6601618158  # twin-prime constant
TWO_FLOOR = 2 * C2 - 1  # ~0.3203  (relative arithmetic floor, §385)

BAND = (14.0, 44.0)
W = 3.0

def lmin_with_infl(cs, w, p_infl, eps):
    """build B_analytic (pole+arch-prime) with Lambda(p_infl)-power weights scaled by eps.
    lam_scale is a MULTIPLIER over NS (build_forms: weight = LAMW * lam_scale). So ones, eps at powers."""
    scale = np.ones_like(L.LAMW)
    for i, n in enumerate(L.NS):
        m = n
        while m % p_infl == 0:
            m //= p_infl
        if m == 1:  # n is a pure power of p_infl
            scale[i] = eps
    f = L.build_forms(cs, w, lam_scale=scale)
    return float(np.linalg.eigvalsh(f['full'])[0])

# pick a well-conditioned sub-edge density: from §411, s=4.0 (n=8, rcond ~6e1, true lmin +8.6e-5, tiny+)
# use s=4.2 where true is +2.2e-3 and clearly positive but near enough the edge to be sensitive.
def find_eps_star(cs, w, p):
    lo, hi = 1.0, 1.30
    # confirm eps=hi breaks it (negative); if not, extend
    if lmin_with_infl(cs, w, p, hi) >= 0:
        hi = 2.0
    # bisect for the sign change
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if lmin_with_infl(cs, w, p, mid) >= 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

print(f"§412 single-prime stiffness. band {BAND} w={W}. two-floor 2C2-1 = {TWO_FLOOR:.4f}")
print(f"{'s':>5} {'n':>3} {'rcond':>8} {'true_lmin':>11}   " + "  ".join(f"eps*-1(p={p})" for p in (2,3,5,7)))
for s in (4.5, 4.2, 4.0, 3.8):
    cs = L.probes(BAND, s)
    if len(cs) < 4: continue
    f0 = L.build_forms(cs, W)
    G = f0['overlap']
    ev = np.linalg.eigvalsh(G)
    rcond = ev[0] / ev[-1]
    tl = float(np.linalg.eigvalsh(f0['full'])[0])
    margins = []
    for p in (2, 3, 5, 7):
        es = find_eps_star(cs, W, p)
        margins.append(es - 1.0)
    print(f"{s:5.1f} {len(cs):3d} {1/rcond:8.1e} {tl:+11.3e}   " +
          "  ".join(f"{m:+.5f}" for m in margins), flush=True)

print(f"""
reading:
  If eps*-1 is roughly PRIME-UNIFORM and ~ {TWO_FLOOR:.3f}, the analytic stiffness = the arithmetic
    two-floor constant 2C2-1 (mechanism made quantitative: the Parity-Bar floors the perturbation).
  If eps*-1 is STRONGEST (smallest) at p=2 and grows with p, the resonance is the low-prime/mean term
    (matches §411 that Lambda(2) is the sensitive direction) -> stiffness is a WEIGHTED single-prime bound
    not a flat constant; the weight Lambda(p)/sqrt(p) is the natural normalizer.
  Either way: a FINITE positive eps* at a well-conditioned s = an unconditional local margin the true
    Lambda enjoys and the second space must reproduce. eps* -> 0 as s -> edge = the margin closes exactly
    at the density where RvM saturates (the sign half lives precisely at the counting edge).
""")

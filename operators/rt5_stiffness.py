#!/usr/bin/env python3
"""§379: ARITHMETIC STIFFNESS — how far can Lambda be bent before the corner margin goes NEG?

§372b: bending Lambda(2) by +5% sends the critical-density margin negative while the TRUE primes
hold it >= 0. Quantify the safety margin: for each band, at the near-critical row, find the
critical multiplicative bend b* on the powers of 2 (Lambda(2^k) -> b*Lambda) at which lmin(full)
first crosses 0. b*-1 = the arithmetic stiffness (how much the true arithmetic can be perturbed
before positivity fails). RH-robustness question: is b*-1 bounded BELOW uniformly across bands
(primes have uniform stiffness), or does it -> 0 (marginal everywhere, no cushion)? Also bend a
NON-arithmetic control (all Lambda by a common factor) to separate "arithmetic-specific" stiffness
from overall scale. Gaussian instrument; the direction of the bend that MOST easily flips it =
the resonant prime (n* by projection).
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def N_rvm(T):
    if T <= 2 * np.pi: return 0.0
    return T / (2 * np.pi) * np.log(T / (2 * np.pi)) - T / (2 * np.pi) + 7.0 / 8.0

def crit_bend(cs, w, which):
    """find bend b on the given prime power's Lambda at which lmin(full) crosses 0.
    which: integer prime power n to bend (all its powers), or 'all' for global scale.
    returns b*-1 (signed smallest-magnitude bend to reach 0); +/- search."""
    base = np.ones_like(L.LAMW)
    def scale_for(b):
        s = np.ones_like(L.LAMW)
        if which == 'all':
            s[:] = b
        else:
            for i, n in enumerate(L.NS):
                m = n
                while m % which == 0: m //= which
                if m == 1: s[i] = b   # powers of `which`
        return s
    def lmin_at(b):
        f = L.build_forms(cs, w, lam_scale=scale_for(b))
        return float(np.linalg.eigvalsh(f['full'])[0])
    l0 = lmin_at(1.0)
    if l0 <= 0: return 0.0, l0
    # increasing b raises the prime term (subtracted) -> lowers margin. search b>1 for zero-crossing.
    lo, hi = 1.0, 1.0
    for _ in range(60):
        hi *= 1.05
        if lmin_at(hi) <= 0: break
    else:
        return None, l0
    for _ in range(50):
        mid = (lo + hi) / 2
        if lmin_at(mid) > 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2 - 1.0, l0

w = 3.0
BANDS = [(14.0, 44.0), (30.0, 60.0), (50.0, 80.0), (75.0, 105.0), (100.0, 130.0)]
print("§379 arithmetic stiffness: critical bend b*-1 to flip corner margin negative (w=3, near-crit row)")
print(f"{'band':>10} {'n':>3} {'lmin0':>10} {'bend n=2':>9} {'bend n=3':>9} {'bend n*':>9} {'bend ALL':>9}")
for band in BANDS:
    b0, b1 = band
    set_grid(b1 + 30.0)
    ncrit = max(4, int(round(N_rvm(b1 + 2 * w) - N_rvm(b0 - 2 * w))) - 1)  # one below crit (margin>0)
    cs = np.linspace(b0, b1, ncrit)
    f = L.build_forms(cs, w)
    # resonant n* by projection
    lA, UA = np.linalg.eigh(f['archpole']); v = UA[:, 0]
    def pform(n):
        u = np.log(n); M = np.zeros((ncrit, ncrit))
        for j in range(ncrit):
            for k in range(j, ncrit):
                M[j, k] = M[k, j] = 2.0 * L.g_pair(cs[j], cs[k], w, u)
        return M
    proj = {}
    for i, n in enumerate(L.NS):
        if L.LAMW[i] == 0 or L.LAMW[i] * np.exp(-(np.log(n)**2)*w*w/4) < 1e-14: continue
        proj[int(n)] = abs(float(v @ pform(int(n)) @ v)) * L.LAMW[i]
    nstar = max(proj, key=lambda n: proj[n]) if proj else 2
    b2, l0 = crit_bend(cs, w, 2)
    b3, _ = crit_bend(cs, w, 3)
    bns, _ = crit_bend(cs, w, nstar)
    ball, _ = crit_bend(cs, w, 'all')
    fmt = lambda x: f"{x:+.4f}" if x is not None else "  none "
    print(f"{str([int(b0),int(b1)]):>10} {ncrit:3d} {l0:10.2e} {fmt(b2)} {fmt(b3)} {fmt(bns)}(n*={nstar}) {fmt(ball)}",
          flush=True)
print("\nreading: b*-1 = arithmetic cushion. If bounded below uniformly (band-flat, >0), primes have")
print("uniform stiffness = quantitative RH-robustness on this front. If ->0 with band, the corner is")
print("marginal everywhere (no cushion, §335 'arithmetic only tested at criticality' made quantitative).")
print("n* being the cheapest to flip = the resonant prime is the RH-critical arithmetic dof.")

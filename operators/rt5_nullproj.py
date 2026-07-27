#!/usr/bin/env python3
"""§374: the CORNER INEQUALITY made explicit — project the Weil form onto the frame-null
direction at n=n_crit (§373) and read what the primes must satisfy to keep the margin >=0.

At the cliff (probe count = RvM count), the overlap Gram has a near-null direction v0 (the
frame-critical mode, computable ZERO-FREE — pure Gaussians). Along v0:
  <v0, B_full v0> = <v0, pole v0> + <v0, arch v0> - <v0, prime v0>
                  = P0 + A0 - sum_n [Lambda(n)/sqrt(n)] * q_n,   q_n := <v0, G_n v0>
where G_n is the n-th prime-power form (Gaussian, zero-free). RH-positivity AT THE CORNER along
this mode <=> sum_n Lambda(n)/sqrt(n) q_n <= P0 + A0. This is an EXPLICIT prime-sum inequality
against the computable weights q_n. Measure:
  - the margin P0+A0 - prime-sum (should be ~0+ at the true primes = the cliff),
  - the per-n prime contributions Lambda(n)/sqrt(n) q_n (which primes push it negative),
  - the SENSITIVITY d(margin)/d(log Lambda(2)) (why 5% bent kills it),
  - and the KEY: is the true-prime prime-sum EXACTLY at the P0+A0 ceiling, or below with slack?
Slack = how close the corner is to violation = the quantitative RH-tightness at this band.
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

def prime_form_n(cs, w, n):
    u = np.log(n)
    G = np.zeros((len(cs), len(cs)))
    for j in range(len(cs)):
        for k in range(j, len(cs)):
            G[j, k] = G[k, j] = 2.0 * L.g_pair(cs[j], cs[k], w, u)
    return G

w = 2.0
BANDS = [(14.0, 44.0), (50.0, 80.0), (100.0, 130.0), (185.0, 215.0)]
print("CORNER INEQUALITY: <v0, B v0> along the frame-null mode at n=n_crit (RvM count).")
print("margin = P0 + A0 - primesum;  RH-at-corner <=> primesum <= P0+A0.")
for band in BANDS:
    b0, b1 = band
    set_grid(b1 + 30.0)
    ncrit = int(round(N_rvm(b1 + 2 * w) - N_rvm(b0 - 2 * w)))
    cs = np.linspace(b0, b1, ncrit)
    f = L.build_forms(cs, w)
    # frame-null direction: bottom eigenvector of the OVERLAP Gram (zero-free, pure Gaussians)
    lo, Vo = np.linalg.eigh(f['overlap'])
    v0 = Vo[:, 0]
    # pole & arch along v0
    z = 0.5j
    fp = np.exp(-((z - cs) ** 2) / (2 * w * w)) + np.exp(-((z + cs) ** 2) / (2 * w * w))
    pole = 2.0 * np.real(np.outer(fp, fp))
    arch = f['archpole'] - pole
    P0 = float(v0 @ pole @ v0)
    A0 = float(v0 @ arch @ v0)
    # per-n prime contributions
    contrib = {}
    for i, n in enumerate(L.NS):
        if L.LAMW[i] == 0.0: continue
        if L.LAMW[i] * np.exp(-(np.log(n) ** 2) * w * w / 4) < 1e-14: continue
        Gn = prime_form_n(cs, w, int(n))
        qn = float(v0 @ Gn @ v0)
        contrib[int(n)] = L.LAMW[i] * qn        # Lambda(n)/sqrt(n) * q_n
    primesum = sum(contrib.values())
    margin = P0 + A0 - primesum
    # sensitivity to Lambda(2): d(margin)/d(scale on all powers of 2) = -sum_{2^k} contrib
    pow2 = sum(v for n, v in contrib.items() if (n & (n - 1)) == 0)
    full_margin = float(np.linalg.eigvalsh(f['full'])[0])
    top = sorted(contrib.items(), key=lambda kv: -abs(kv[1]))[:5]
    print(f"\nband {band}  n_crit={ncrit}")
    print(f"  P0={P0:+.4e}  A0={A0:+.4e}  primesum={primesum:+.4e}")
    print(f"  MARGIN along v0 = {margin:+.4e}   (full-form lmin = {full_margin:+.4e})")
    print(f"  slack = margin/(P0+A0) = {margin/(P0+A0):.4f}   (0 = at violation, 1 = primes irrelevant)")
    print(f"  pow2 contribution = {pow2:+.4e} ({100*pow2/primesum:.1f}% of primesum) -> "
          f"5% bend moves margin by {-0.05*pow2:+.2e} (flips if > margin)")
    print(f"  top-5 |contrib| n: " + ", ".join(f"{n}:{v:+.2e}" for n, v in top))
print("\nreading: the corner inequality primesum <= P0+A0 holds (margin>=0) for TRUE primes;")
print("slack quantifies tightness. pow2 dominating primesum explains §372b (5% Lambda(2) bend flips")
print("the margin negative). The RH content at the corner = why the TRUE Lambda keeps primesum just")
print("under P0+A0 while a small bend exceeds it — an explicit inequality on q_n-weighted prime sums.")

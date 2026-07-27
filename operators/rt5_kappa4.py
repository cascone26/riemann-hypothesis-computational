#!/usr/bin/env python3
"""§384: instrument (a) — connected kappa_4 / excess-kurtosis gamma_2 of the corner statistic
Sigma c_w delta-Lambda_w, real vs PAIR-MATCHED (phase-randomized) surrogate. Pre-registered §383.

Sigma c_w dL_w is linear in the prime fluctuations dL = Lambda-1; its gamma_2 = kappa_4/kappa_2^2
vanishes for ANY Gaussian/covariance (pair) model => gamma_2 != 0 is a STRUCTURAL beyond-pair
certificate. Test: block the log-domain prime fluctuation up the line, corner-weight it (Gaussian
bump at the resonant scale = the §374 corner mode's spectral content), measure gamma_2 across
blocks. Pair-matched null = PHASE-RANDOMIZED surrogate (preserves 2-point/power spectrum exactly,
Gaussianizes higher cumulants) — the textbook pair-matched ensemble. 3-way falsifier on sign(gamma_2):
  gamma_2 < 0 (clearing 0, |g2| covers O(1)) => platykurtic => Sigma c_w >= 0 pinned = 3-clause object
  gamma_2 > 0 => beyond-pair CERTIFIED, RH-pin disproved
  gamma_2 ~ 0 => pair-dominated at 4th order, wall re-lettered
HONEST LABEL: this measures gamma_2 of the corner-weighted prime-fluctuation sum; the corner weight
is a log-Gaussian proxy for the §374 q_n profile. Surrogate = phase randomization (standard).
"""
import numpy as np

rng = np.random.default_rng(20260728)
NMAX = 2_000_000
print(f"sieving Lambda to {NMAX} ...", flush=True)
lam = np.zeros(NMAX + 1)
comp = np.zeros(NMAX + 1, dtype=bool)
for p in range(2, NMAX + 1):
    if not comp[p]:
        comp[2 * p::p] = True
        lp = np.log(p); pk = p
        while pk <= NMAX:
            lam[pk] = lp; pk *= p

# log-domain fluctuation field on a uniform u=log n grid
n = np.arange(2, NMAX + 1)
u = np.log(n)
dL = lam[2:] - 1.0                        # Lambda(n) - 1  (mean-zero by PNT)
# resample onto a uniform u-grid (bin-average) so phase randomization is well-defined
UMAX = np.log(NMAX)
NG = 1 << 19
ubins = np.linspace(np.log(2), UMAX, NG + 1)
idx = np.clip(np.searchsorted(ubins, u) - 1, 0, NG - 1)
field = np.zeros(NG); cnt = np.zeros(NG)
np.add.at(field, idx, dL); np.add.at(cnt, idx, 1.0)
field = np.where(cnt > 0, field / np.maximum(cnt, 1), 0.0)
field = field - field.mean()

def excess_kurt(x):
    x = x - x.mean(); m2 = np.mean(x * x); m4 = np.mean(x ** 4)
    return m4 / m2 ** 2 - 3.0

def phase_randomize(x):
    X = np.fft.rfft(x)
    ph = np.exp(1j * rng.uniform(0, 2 * np.pi, len(X)))
    ph[0] = 1.0
    if len(x) % 2 == 0: ph[-1] = 1.0
    return np.fft.irfft(X * np.abs(1) * ph, n=len(x))   # preserve |X|, randomize phase

du = ubins[1] - ubins[0]
print(f"{'u*':>5} {'sigma':>6} {'g2_real':>9} {'g2_surr mean':>13} {'g2_surr sd':>11} {'z':>7} {'verdict':>10}")
# corner weight = log-Gaussian bump at resonant scale u* (Bragg ~ 2pi/s; scan a few), width sigma
for u_star, sigma in [(3.0, 0.8), (4.0, 1.0), (5.0, 1.2), (6.0, 1.5)]:
    w = np.exp(-((np.arange(NG) * du + np.log(2)) - u_star) ** 2 / (2 * sigma ** 2))
    g = field * w                         # corner-weighted fluctuation
    # block into windows; block-sum = Sigma c_w dL_w per window
    BL = 2048
    nb = NG // BL
    blocks = g[:nb * BL].reshape(nb, BL).sum(axis=1)
    g2_real = excess_kurt(blocks)
    # phase-randomized surrogates of the weighted field, re-blocked
    surr = []
    for _ in range(40):
        gs = phase_randomize(g)
        bs = gs[:nb * BL].reshape(nb, BL).sum(axis=1)
        surr.append(excess_kurt(bs))
    surr = np.array(surr)
    z = (g2_real - surr.mean()) / (surr.std() + 1e-12)
    verd = ("BEYOND-PAIR+" if z > 3 else "BEYOND-PAIR-" if z < -3 else "pair~0")
    print(f"{u_star:5.1f} {sigma:6.2f} {g2_real:+9.3f} {surr.mean():+13.3f} {surr.std():11.3f} "
          f"{z:+7.2f} {verd:>12}", flush=True)
print("\nreading: sign(g2_real - surr) with |z|>3 = structural beyond-pair (surrogate kills all pair");
print("structure, so any robust excess is connected kappa_4). NEGATIVE => platykurtic => the pin lead;")
print("POSITIVE => beyond-pair real but pin disproved; ~0 => pair-dominated at 4th order.")

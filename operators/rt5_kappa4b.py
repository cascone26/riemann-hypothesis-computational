#!/usr/bin/env python3
"""§384b: instrument (a) DONE RIGHT — IAAFT pair-matched surrogate (preserves marginal AND power
spectrum, Gaussianizes only higher correlations). gamma_2_real - gamma_2_IAAFT isolates the
CONNECTED kappa_4 (beyond-pair correlation), immune to the §384 marginal-spikiness confound.
3-way falsifier on sign (pre-reg §383): <0 platykurtic pin lead; >0 beyond-pair, pin disproved;
~0 pair-dominated / wall re-lettered."""
import numpy as np
rng = np.random.default_rng(20260728)
NMAX = 2_000_000
print(f"sieving Lambda to {NMAX} ...", flush=True)
lam = np.zeros(NMAX + 1); comp = np.zeros(NMAX + 1, dtype=bool)
for p in range(2, NMAX + 1):
    if not comp[p]:
        comp[2 * p::p] = True
        lp = np.log(p); pk = p
        while pk <= NMAX: lam[pk] = lp; pk *= p
n = np.arange(2, NMAX + 1); u = np.log(n); dL = lam[2:] - 1.0
UMAX = np.log(NMAX); NG = 1 << 19
ubins = np.linspace(np.log(2), UMAX, NG + 1)
idx = np.clip(np.searchsorted(ubins, u) - 1, 0, NG - 1)
field = np.zeros(NG); cnt = np.zeros(NG)
np.add.at(field, idx, dL); np.add.at(cnt, idx, 1.0)
field = np.where(cnt > 0, field / np.maximum(cnt, 1), 0.0); field -= field.mean()
du = ubins[1] - ubins[0]

def excess_kurt(x):
    x = x - x.mean(); m2 = np.mean(x * x); return np.mean(x ** 4) / m2 ** 2 - 3.0

def iaaft(x, iters=100):
    # Schreiber-Schmitz: preserve amplitude spectrum |X| AND sorted marginal
    xs = np.sort(x); amp = np.abs(np.fft.rfft(x))
    y = rng.permutation(x)
    for _ in range(iters):
        Y = np.fft.rfft(y); Y = amp * np.exp(1j * np.angle(Y))
        y = np.fft.irfft(Y, n=len(x))
        r = np.argsort(np.argsort(y)); y = xs[r]     # impose marginal
    return y

print(f"{'u*':>5} {'sig':>5} {'g2_real':>9} {'g2_iaaft mean':>13} {'sd':>7} {'z':>7} {'verdict':>13}")
for u_star, sigma in [(3.0, 0.8), (4.0, 1.0), (5.0, 1.2)]:
    w = np.exp(-((np.arange(NG) * du + np.log(2)) - u_star) ** 2 / (2 * sigma ** 2))
    g = field * w
    BL = 2048; nb = NG // BL
    g2_real = excess_kurt(g[:nb * BL].reshape(nb, BL).sum(1))
    surr = []
    for _ in range(20):
        gs = iaaft(g)
        surr.append(excess_kurt(gs[:nb * BL].reshape(nb, BL).sum(1)))
    surr = np.array(surr); z = (g2_real - surr.mean()) / (surr.std() + 1e-12)
    verd = ("BEYOND-PAIR+ (pin disproved)" if z > 3 else
            "BEYOND-PAIR- (PIN LEAD)" if z < -3 else "pair-dominated ~0")
    print(f"{u_star:5.1f} {sigma:5.2f} {g2_real:+9.3f} {surr.mean():+13.3f} {surr.std():7.3f} {z:+7.2f} {verd:>13}",
          flush=True)
print("\nreading: IAAFT preserves marginal+spectrum => g2_real-g2_iaaft is the CONNECTED kappa_4 only");
print("(marginal spikiness now in BOTH, cancels). sign+|z|>3 = clean beyond-pair verdict per §383.")

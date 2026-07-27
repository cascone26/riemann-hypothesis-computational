#!/usr/bin/env python3
"""§384c: kappa_4 resolver — larger N (8e6), high scales u*=5,6,7, resolve the §384b z-trend.
z stable/up with N => real emerging beyond-pair; z down => few-sample noise at sparse large-u."""
import numpy as np
rng = np.random.default_rng(7)
NMAX = 8_000_000
print(f"sieving to {NMAX} ...", flush=True)
lam = np.zeros(NMAX + 1); comp = np.zeros(NMAX + 1, dtype=bool)
for p in range(2, NMAX + 1):
    if not comp[p]:
        comp[2 * p::p] = True; lp = np.log(p); pk = p
        while pk <= NMAX: lam[pk] = lp; pk *= p
n = np.arange(2, NMAX + 1); u = np.log(n); dL = lam[2:] - 1.0
NG = 1 << 20; ub = np.linspace(np.log(2), np.log(NMAX), NG + 1)
idx = np.clip(np.searchsorted(ub, u) - 1, 0, NG - 1)
field = np.zeros(NG); cnt = np.zeros(NG)
np.add.at(field, idx, dL); np.add.at(cnt, idx, 1.0)
field = np.where(cnt > 0, field / np.maximum(cnt, 1), 0.0); field -= field.mean()
du = ub[1] - ub[0]
def ek(x): x = x - x.mean(); m2 = np.mean(x * x); return np.mean(x ** 4) / m2 ** 2 - 3
def iaaft(x, it=100):
    xs = np.sort(x); amp = np.abs(np.fft.rfft(x)); y = rng.permutation(x)
    for _ in range(it):
        Y = np.fft.rfft(y); Y = amp * np.exp(1j * np.angle(Y)); y = np.fft.irfft(Y, n=len(x))
        y = xs[np.argsort(np.argsort(y))]
    return y
print(f"{'u*':>4} {'g2_real':>9} {'iaaft':>8} {'sd':>6} {'z':>7} {'nblk':>6}", flush=True)
for u_star, sig in [(5.0, 1.2), (6.0, 1.4), (7.0, 1.6)]:
    w = np.exp(-((np.arange(NG) * du + np.log(2)) - u_star) ** 2 / (2 * sig ** 2))
    g = field * w; BL = 1024; nb = NG // BL
    gr = ek(g[:nb * BL].reshape(nb, BL).sum(1))
    s = np.array([ek(iaaft(g)[:nb * BL].reshape(nb, BL).sum(1)) for _ in range(40)])
    z = (gr - s.mean()) / (s.std() + 1e-12)
    print(f"{u_star:4.1f} {gr:+9.3f} {s.mean():+8.3f} {s.std():6.3f} {z:+7.2f} {nb:6d}", flush=True)
print("resolver done: vs §384b(NMAX=2e6) u*=5 z=+2.79 — up=>real, down=>few-sample noise.", flush=True)

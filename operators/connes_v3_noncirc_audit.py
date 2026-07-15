#!/usr/bin/env python3
"""
Non-circularity audit for connes_decisive_v3.

The v3 build never sees a single zeta zero (only digamma Omega, von Mangoldt,
W_{0,2}).  Decisive test of genuineness: dump EVERY root of xi_hat in a wide
window and compare against INDEPENDENT ground truth from mpmath.zetazero --
including high zeros (#13..#20) that appear NOWHERE in any source file.  If the
model emits gamma_13 ~ 59.347 etc. that it was never given, the reproduction is
a genuine non-circular prediction, not a planted constant or a 2pi/L lattice.
"""
import numpy as np, mpmath as mp
from connes_decisive_v3 import build_tau, xi_hat_factory, find_zeros, L, TWO_PI

mp.mp.dps = 25
# independent ground truth, computed here, never fed into the construction
GT = [float(mp.im(mp.zetazero(k))) for k in range(1, 21)]

for N in [32, 48, 64]:
    R = TWO_PI * N / L + 80.0
    G = max(60000, 3000 * N)
    tau, idx = build_tau(N, R, G)
    ev, evec = np.linalg.eigh(tau)
    k0 = int(np.argmin(np.abs(ev)))
    v = evec[:, k0]
    D = 2 * N + 1
    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])
    s = even.sum() / np.sqrt(L)
    if abs(s) > 1e-30:
        even = even / s
    xh = xi_hat_factory(even, idx)
    zeros = find_zeros(xh, 8, 75, 6000)
    print(f"\n=== N={N}  eps_min={ev[0]:.3e} eps_kernel={ev[k0]:.3e}  "
          f"{len(zeros)} roots in [8,75] ===")
    # lattice check: 2pi/L spacing would be ~2.3808
    if len(zeros) > 1:
        diffs = np.diff(sorted(zeros))
        print(f"  root spacing  min={diffs.min():.4f} "
              f"max={diffs.max():.4f}  (2pi/L lattice = {TWO_PI/L:.4f})")
    for k in range(20):
        g = GT[k]
        if zeros:
            got = min(zeros, key=lambda z: abs(z - g))
            err = abs(got - g)
            seen = "in TRUE_ZEROS" if k < 12 else "NOT in any file"
            fl = "OK" if err < 1e-3 else "  "
            print(f"  gamma_{k+1:<2d} {g:10.6f}  model {got:10.6f}  "
                  f"err {err:.2e} {fl}   [{seen}]")

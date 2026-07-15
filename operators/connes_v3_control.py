#!/usr/bin/env python3
"""
Control / sensitivity probe for connes_decisive_v3.

The archimedean symbol Omega(r)=Re psi(1/4 + i r/2) - log pi has the '1/4'
ENFORCED by the critical line: it is s/2 at s = 1/2 + i r.  If the construction
genuinely encodes RH positivity, replacing 1/4 by a != 1/4 (i.e. pretending the
zeros sit on Re s = 2a) MUST break it: eps_N goes strongly negative and the
secular zeros stop being the zeta zeros.  If instead any 'a' still yields a PSD
tau and the same zeros, the v3 signal is a normalization artifact and vacuous.

Only the structural symbol changes; basis, primes, boundary are untouched.
"""
import numpy as np, mpmath as mp
from scipy.special import digamma as _digamma
import connes_decisive_v3 as v3

mp.mp.dps = 25
GT = [float(mp.im(mp.zetazero(k))) for k in range(1, 7)]   # independent truth


def build_tau_a(N, R, G, a):
    idx = list(range(-N, N + 1))
    D = len(idx)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = np.real(_digamma(a + 0.5j * grid)) - v3.LOGPI      # only change: a vs 1/4
    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = v3.e_hat_vec(grid, n)
    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)
    tau = arch + v3.w02_mat(idx) - v3.prime_mat(idx)
    return (tau + tau.T) / 2, idx


N = 32
R = v3.TWO_PI * N / v3.L + 80.0
G = 3000 * N
print(f"N={N}  (a=0.25 is the true critical-line symbol s/2 at s=1/2+ir)\n")
for a in [0.10, 0.25, 0.40, 0.60]:
    tau, idx = build_tau_a(N, R, G, a)
    ev, evec = np.linalg.eigh(tau)
    k0 = int(np.argmin(np.abs(ev)))
    v = evec[:, k0]
    D = 2 * N + 1
    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])
    s = even.sum() / np.sqrt(v3.L)
    if abs(s) > 1e-30:
        even = even / s
    xh = v3.xi_hat_factory(even, idx)
    zeros = v3.find_zeros(xh, 8, 45, 4000)
    errs = []
    for g in GT:
        if zeros:
            errs.append(abs(min(zeros, key=lambda z: abs(z - g)) - g))
        else:
            errs.append(float("nan"))
    maxe = max(errs)
    tag = "  <-- TRUE critical-line symbol" if abs(a - 0.25) < 1e-9 else ""
    verdict = "ZEROS REPRODUCED" if maxe < 1e-3 else "broken (zeros scatter)"
    print(f"a={a:.2f}  eps_min={ev[0]:+.3e}  eps_kernel={ev[k0]:+.3e}  "
          f"max|err γ1..6|={maxe:.2e}  {verdict}{tag}")

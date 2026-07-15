#!/usr/bin/env python3
"""
connes_decisive_v3.py -- v2's NON-circular construction, vectorized so it runs.

Same math as v2 (correct, non-circular):
   tau_{n,m} = (1/2pi) int Omega(r) e_hat_n(r) conj(e_hat_m(r)) dr   [arch]
             + W_{0,2}(n,m)                                          [boundary/pole]
             - sum_{k<=lambda^2} Lambda(k) k^{-1/2} q(U_n,U_m)(log k) [prime]
with Omega(r) = Re psi(1/4 + i r/2) - log pi   (digamma; ZERO-free; no zeta zeros in).

v2 was correct but called mpmath.digamma (dps=30) ~10^7 times -> never finished.
v3: ONE shared dense grid per N; Omega via scipy.special.digamma (vectorized,
complex); the whole archimedean block is a single matrix product
   arch = (1/2pi) * E @ diag(w * Omega) @ E^H
E[k,g] = e_hat_{idx[k]}(grid[g]), w = trapezoid weights.  Identical integrand,
~10000x faster, no precision loss that matters for a float Gram matrix.

Decisive test: eps_N -> 0^+ (Weil positivity; v1 gave -3.6) AND zeros of xi_hat
snap to the true zeta zeros instead of the 2pi/L Dirichlet lattice.
"""

import json, time
import numpy as np
import mpmath as mp
from scipy.special import digamma as _digamma

mp.mp.dps = 30
LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))

TRUE_ZEROS = [14.134725142, 21.022039639, 25.010857580, 30.424876126,
              32.935061588, 37.586178159, 40.918719012, 43.327073281,
              48.005150881, 49.773832478, 52.970321478, 56.446247697]


def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab
VM = vm_table(LAMBDA_SQ)


def omega_vec(r):
    """Omega(r) = Re psi(1/4 + i r/2) - log pi, vectorized over array r."""
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI


def e_hat_vec(r, n):
    """e_hat_n(r) = (1/sqrt L) e^{i d L/2} 2 sin(d L/2)/d, d = s_n - r, vectorized."""
    s = TWO_PI * n / L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L) + 0j
    dd = d[~small]
    out[~small] = (1.0 / np.sqrt(L)) * np.exp(0.5j * dd * L) * 2.0 * np.sin(0.5 * dd * L) / dd
    return out


def q_kernel(n, m, y):
    if n == m:
        return 2.0 * (1.0 - abs(y) / L) * np.cos(TWO_PI * n * y / L)
    return (np.sin(TWO_PI * m * y / L) - np.sin(TWO_PI * n * y / L)) / (np.pi * (n - m))


def prime_mat(idx):
    D = len(idx)
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            s = 0.0
            for k, Lam in VM:
                s += Lam * k ** (-0.5) * q_kernel(n, m, np.log(k))
            P[ii, jj] = P[jj, ii] = s
    return P


def w02_mat(idx):
    D = len(idx)
    W = np.zeros((D, D))
    sh = np.sinh(L / 4) ** 2
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = (L ** 2 - 16 * np.pi ** 2 * m * n)
            den = (L ** 2 + 16 * np.pi ** 2 * m ** 2) * (L ** 2 + 16 * np.pi ** 2 * n ** 2)
            W[ii, jj] = 32 * L * sh * num / den
    return W


def build_tau(N, R, G):
    idx = list(range(-N, N + 1))
    D = len(idx)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5 * h          # trapezoid weights
    Om = omega_vec(grid)
    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = e_hat_vec(grid, n)
    A = E * (w * Om)                                    # (D,G)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)        # (D,D)
    tau = arch + w02_mat(idx) - prime_mat(idx)
    tau = (tau + tau.T) / 2
    return tau, idx


def xi_hat_factory(xi, idx):
    poles = np.array([TWO_PI * j / L for j in idx])
    xi = np.array(xi, dtype=float)
    def xh(z):
        d = z - poles
        d[np.abs(d) < 1e-12] = 1e-12
        return (2.0 / np.sqrt(L)) * np.sin(z * L / 2.0) * np.sum(xi / d)
    return xh


def find_zeros(xh, lo, hi, ng):
    xs = np.linspace(lo, hi, ng)
    vs = np.array([xh(x) for x in xs])
    out = []
    for i in range(1, len(xs)):
        if vs[i - 1] == 0 or (vs[i - 1] < 0) != (vs[i] < 0):
            a, b = xs[i - 1], xs[i]
            for _ in range(80):
                mid = 0.5 * (a + b)
                if (xh(a) < 0) != (xh(mid) < 0): b = mid
                else: a = mid
            r = 0.5 * (a + b)
            if not out or abs(r - out[-1]) > 1e-4:
                out.append(r)
    return out


def score(zeros):
    res = []
    for tz in TRUE_ZEROS:
        if zeros:
            g = min(zeros, key=lambda z: abs(z - tz))
            res.append((tz, g, abs(g - tz)))
        else:
            res.append((tz, None, None))
    return res


def run(N):
    R = TWO_PI * N / L + 80.0
    G = max(60000, 3000 * N)
    t0 = time.time()
    tau, idx = build_tau(N, R, G)
    ev, evec = np.linalg.eigh(tau)
    k0 = int(np.argmin(np.abs(ev)))
    eps_min = float(ev[0]); eps_kernel = float(ev[k0])
    v = evec[:, k0]
    D = 2 * N + 1
    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])
    s = even.sum() / np.sqrt(L)
    if abs(s) > 1e-30: even = even / s
    xh = xi_hat_factory(even, idx)
    zeros = find_zeros(xh, 10, 60, 3000)
    return eps_min, eps_kernel, zeros, time.time() - t0


if __name__ == "__main__":
    out = {"lambda_sq": LAMBDA_SQ, "L": L,
           "construction": "v3 vectorized: arch = (1/2pi) E diag(w Omega) E^H, "
           "Omega=Re psi(1/4+ir/2)-log pi (scipy digamma, zero-free); "
           "prime = von Mangoldt; boundary = W_{0,2}. No zeta zeros as input.",
           "runs": []}
    for N in [16, 24, 32, 48]:
        eps_min, eps_k, zeros, dt = run(N)
        sc = score(zeros)
        print(f"\n=== N={N}  eps_min={eps_min:.4e}  eps_kernel={eps_k:.4e}  "
              f"({dt:.1f}s)  {len(zeros)} zeros[10,60] ===")
        for tz, g, e in sc:
            gs = f"{g:.5f}" if g is not None else "  --  "
            es = f"{e:.2e}" if e is not None else " n/a "
            fl = "OK" if (e is not None and e < 0.05) else "  "
            print(f"  true {tz:9.5f}  got {gs}  err {es} {fl}")
        out["runs"].append({"N": N, "eps_min": eps_min, "eps_kernel": eps_k,
                            "secs": dt, "n_zeros": len(zeros), "zeros": zeros,
                            "score": [(a, b, c) for a, b, c in sc]})
    p = "/Users/scones/projects/riemann/results/connes_decisive_v3.json"
    json.dump(out, open(p, "w"), indent=2)
    print(f"\nsaved -> {p}")

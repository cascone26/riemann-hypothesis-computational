#!/usr/bin/env python3
"""
connes_decisive_v2.py -- the NON-CIRCULAR archimedean term, done correctly.

v1 mistake: reused the prime q-kernel (sin-difference / (1-|y|/L)cos) for the
archimedean term and paired it with a guessed Hadamard-FP density.  The doc
(sec 4.1) says W_R needs x*sin, x*cos terms too -> that is the digamma symbol,
not rho-against-q.  Result: tau indefinite, eps_N -> -3.6, xi collapsed to the
2 pi / L Dirichlet lattice.

THE INSIGHT.  connes_exact_weil.py builds
      tau_arch = sum over zeta zeros gamma of  e_hat_n(gamma) e_hat_m(gamma)*
which is circular.  But the Weil explicit formula says EXACTLY:
   sum_gamma F(gamma)  =  (1/2pi) int Omega(r) F(r) dr  +  pole  -  sum_primes
with the archimedean symbol
   Omega(r) = Re psi(1/4 + i r/2) - log pi          (psi = digamma; NO zeros)
So the correct, non-circular archimedean matrix element is the SAME sinc kernel
e_hat integrated against Omega(r) dr -- zeros never enter.

      tau_{n,m} = (1/2pi) int Omega(r) e_hat_n(r) conj(e_hat_m(r)) dr   [arch]
                +  W_{0,2}(n,m)                                          [boundary/pole]
                -  sum_{k<=lambda^2} Lambda(k) k^{-1/2} q(U_n,U_m)(log k) [prime]

Test: does eps_N now -> 0^+ (Weil positivity) and do the zeros of xi_hat snap
to the true zeta zeros instead of the lattice?
"""

import sys, json, time
import numpy as np
import mpmath as mp

mp.mp.dps = 30

LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))     # = log 14
TWO_PI = 2.0 * np.pi

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

# archimedean symbol  Omega(r) = Re psi(1/4 + i r/2) - log pi   (digamma; zero-free)
LOGPI = float(mp.log(mp.pi))
def Omega(r):
    return float(mp.re(mp.digamma(mp.mpf('0.25') + 0.5j * mp.mpf(r)))) - LOGPI

# sinc kernel  e_hat_n(r) = (1/sqrt L) e^{i (s_n - r) L/2} 2 sin((s_n-r)L/2)/(s_n-r)
def e_hat(r, n):
    s = TWO_PI * n / L
    d = s - r
    if abs(d) < 1e-12:
        return np.sqrt(L) + 0j                 # limit value
    return (1.0 / np.sqrt(L)) * np.exp(0.5j * d * L) * 2.0 * np.sin(0.5 * d * L) / d

# archimedean matrix element: (1/2pi) int Omega(r) e_hat_n(r) conj(e_hat_m(r)) dr
def arch_entry(n, m, R, nq):
    sn, sm = TWO_PI * n / L, TWO_PI * m / L
    # break the real line at the two sinc peaks; dense grid, trapezoid in float
    lo, hi = -R, R
    pts = np.unique(np.concatenate([
        np.linspace(lo, hi, nq),
        np.linspace(min(sn, sm) - 6, max(sn, sm) + 6, 4000)]))
    pts.sort()
    om = np.array([Omega(r) for r in pts])
    en = np.array([e_hat(r, n) for r in pts])
    em = np.array([e_hat(r, m) for r in pts])
    integ = om * en * np.conj(em)
    _trap = getattr(np, "trapezoid", None) or np.trapz
    val = _trap(integ, pts) / (2.0 * np.pi)
    return float(np.real(val))

# prime kernel q(U_n,U_m)(y)
def q_kernel(n, m, y):
    if n == m:
        return 2.0 * (1.0 - abs(y) / L) * np.cos(TWO_PI * n * y / L)
    return (np.sin(TWO_PI * m * y / L) - np.sin(TWO_PI * n * y / L)) / (np.pi * (n - m))

def prime_entry(n, m):
    s = 0.0
    for k, Lam in VM:
        y = np.log(k)
        s += Lam * k ** (-0.5) * q_kernel(n, m, y)
    return s

def W02(n, m):
    num = (L ** 2 - 16 * np.pi ** 2 * m * n)
    den = (L ** 2 + 16 * np.pi ** 2 * m ** 2) * (L ** 2 + 16 * np.pi ** 2 * n ** 2)
    return 32 * L * np.sinh(L / 4) ** 2 * num / den

def build_tau(N, R, nq):
    idx = list(range(-N, N + 1))
    D = 2 * N + 1
    tau = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            v = arch_entry(n, m, R, nq) + W02(n, m) - prime_entry(n, m)
            tau[ii, jj] = v
            tau[jj, ii] = v
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
                fm = xh(mid)
                if (xh(a) < 0) != (fm < 0): b = mid
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
    R = TWO_PI * N / L + 60.0
    nq = max(8000, 400 * N)
    t0 = time.time()
    tau, idx = build_tau(N, R, nq)
    tau = (tau + tau.T) / 2
    ev, evec = np.linalg.eigh(tau)
    k0 = int(np.argmin(np.abs(ev)))          # vector nearest the KERNEL (eps_N ~ 0)
    eps_min = float(ev[0]); eps_kernel = float(ev[k0])
    v = evec[:, k0]
    D = 2 * N + 1
    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])
    s = even.sum() / np.sqrt(L)
    if abs(s) > 1e-30: even = even / s
    xh = xi_hat_factory(even, idx)
    zeros = find_zeros(xh, 10, 60, 2400)
    return eps_min, eps_kernel, zeros, time.time() - t0

if __name__ == "__main__":
    out = {"lambda_sq": LAMBDA_SQ, "L": L,
           "construction": "NON-circular: archimedean = (1/2pi) int Omega(r) "
           "e_hat_n e_hat_m dr, Omega=Re psi(1/4+ir/2)-log pi (digamma, zero-free); "
           "prime = von Mangoldt; boundary = W_{0,2}. No zeta zeros as input.",
           "runs": []}
    for N in [16, 24, 32]:
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
    p = "/Users/scones/projects/riemann/results/connes_decisive_v2.json"
    json.dump(out, open(p, "w"), indent=2)
    print(f"\nsaved -> {p}")

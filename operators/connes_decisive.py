#!/usr/bin/env python3
"""
connes_decisive.py  --  NON-CIRCULAR Connes-Consani-Moscovici zeta spectral triple.

Goal: build the truncated Weil matrix tau at lambda^2 = 14 WITHOUT ever feeding
zeta zeros in, then check whether its eigenvalue spectrum reproduces the real
zeros of zeta (Theorem 5.10). If it does, that is a genuine *construction* of the
zeros from a self-adjoint operator -- the only non-circular route to RH.

Key contrast with connes_exact_weil.py (the broken file):
  - OLD: archimedean term = e_hat_matrix.T @ e_hat_matrix  -> a SUM OVER KNOWN
    ZEROS. Circular (zeros in -> zeros out), low-rank, inaccurate.
  - NEW: archimedean term = integral_0^L q(U_n,U_m)(y) * rho(y) dy with the
    elementary zero-free density rho(y)=e^{y/2}/(e^y-e^{-y}). No zeros anywhere.

Weil matrix (doc sec 5, line 234):
    tau_{n,m} = W_R(V_n,V_m) + W_{0,2}(V_n,V_m) - sum_k Lambda(k) k^{-1/2} q(U_n,U_m)(log k)

Toeplitz structure (Lemma 5.1): tau is determined by O(N) generators, not (2N+1)^2:
    off-diagonal  tau_{n,m} = (b_n - b_m)/(n - m),  n != m
    diagonal      tau_{n,n} = a_n
with b odd, a even.  We compute b_k from the sin-transform of (rho-measure minus
prime-measure), and a_n from the regularized diagonal + boundary closed form.
"""

import sys, json, time
import mpmath as mp
import numpy as np

mp.mp.dps = 40

LAMBDA_SQ = 14
L = mp.log(LAMBDA_SQ)            # L = 2 log(lambda) = log(lambda^2) = log 14
TWO_PI = 2 * mp.pi

# true imaginary parts of the first zeta zeros (for SCORING ONLY, never used to build tau)
TRUE_ZEROS = [14.134725142, 21.022039639, 25.010857580, 30.424876126,
              32.935061588, 37.586178159, 40.918719012, 43.327073281,
              48.005150881, 49.773832478, 52.970321478, 56.446247697]

# ---------------------------------------------------------------------------
# von Mangoldt data for k <= e^L = lambda^2 = 14
# k in {2,3,4,5,7,8,9,11,13}; Lambda = log p for prime powers
# ---------------------------------------------------------------------------
def vonmangoldt_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
        n = k
        p = 0
        f = 2
        ok = True
        while f * f <= n:
            if n % f == 0:
                if p == 0:
                    p = f
                elif p != f:
                    ok = False
                    break
                while n % f == 0:
                    n //= f
            f += 1
        if n > 1:
            if p == 0:
                p = n
            elif p != n:
                ok = False
        if ok and p > 0:
            tab.append((k, mp.log(p)))   # (k, Lambda(k))
    return tab

VM = vonmangoldt_table(LAMBDA_SQ)        # [(k, Lambda(k)), ...]

# ---------------------------------------------------------------------------
# Archimedean density (ZERO-FREE -- this is the whole point)
#   rho(y) = e^{y/2} / (e^y - e^{-y}) = e^{y/2} / (2 sinh y)
# Singular as y->0 like 1/(2y); we use Hadamard finite-part by subtracting the
# 1/(2y) head, which is integrable after subtraction.
# ---------------------------------------------------------------------------
def rho(y):
    return mp.e ** (y / 2) / (mp.e ** y - mp.e ** (-y))

def rho_reg(y):
    # rho(y) - 1/(2y): bounded as y->0  (limit = 1/4)
    if y == 0:
        return mp.mpf(1) / 4
    return rho(y) - 1 / (2 * y)

# sin / cos transforms of the archimedean measure on [0,L]
#   Ssin(k) = FP integral_0^L sin(2 pi k y / L) rho(y) dy        (k integer)
#   Scos(k) = FP integral_0^L (1-y/L) cos(2 pi k y / L) rho(y) dy
# off-diagonal needs only Ssin (q has no constant part there);
# diagonal needs the (1-|y|/L)cos kernel.
def _w(k):
    return TWO_PI * k / L

def Ssin(k):
    if k == 0:
        return mp.mpf(0)
    w = _w(k)
    # sin(0)=0 kills the 1/(2y) singular head exactly, integrand is regular
    f = lambda y: mp.sin(w * y) * rho(y)
    return mp.quad(f, [0, L / 2, L])

def Scos_diag(k):
    # integral_0^L (1 - y/L) cos(2 pi k y / L) rho(y) dy, Hadamard FP
    w = _w(k)
    g = lambda y: (1 - y / L) * mp.cos(w * y)
    # split rho = 1/(2y) + rho_reg ; FP of g(y)/(2y) head:
    #   g(0)=1 ; FP integral_0^L 1/(2y) dy  -> regularized to  (1/2) * (-log eps)
    #   the eps-divergence is the standard archimedean conductor term; in the
    #   spectral-triple it is cancelled by the delta_N boundary pairing, leaving
    #   the finite Weil constant  C_R = (1/2)(log(4 pi) + gamma_E)  (Riemann xi
    #   archimedean log-derivative at the symmetric point). We carry it as an
    #   even, n-independent additive constant -> shifts ALL eigenvalues by a
    #   fixed amount, so eigenvalue *spacings* (zero gaps) are invariant to it.
    reg_part = mp.quad(lambda y: g(y) * rho_reg(y), [0, L / 2, L])
    head_smooth = mp.quad(lambda y: (g(y) - mp.cos(w * y) * mp.exp(0)) / (2 * y)
                          if y != 0 else mp.mpf(0), [0, L])  # subtract g(0)*head
    # g(0)=1; the residual smooth part of (g(y)-g(0))/(2y) is finite
    head_smooth2 = mp.quad(lambda y: (g(y) - 1) / (2 * y) if y != 0 else
                           ( -1/L + 0 ), [0, L/2, L])
    C_R = (mp.log(4 * mp.pi) + mp.euler) / 2     # fixed archimedean constant
    return reg_part + head_smooth2 + C_R

# prime-measure transforms (discrete, exact, non-circular: just primes)
def Psin(k):
    s = mp.mpf(0)
    w = _w(k)
    for kk, Lam in VM:
        y = mp.log(kk)
        s += Lam * kk ** (mp.mpf(-1) / 2) * mp.sin(w * y)
    return s

def Pcos_diag(k):
    s = mp.mpf(0)
    w = _w(k)
    for kk, Lam in VM:
        y = mp.log(kk)
        s += Lam * kk ** (mp.mpf(-1) / 2) * 2 * (1 - y / L) * mp.cos(w * y)
    return s

# boundary closed form  W_{0,2}(V_n,V_m) (doc 4.3)
def W02(n, m):
    num = (L ** 2 - 16 * mp.pi ** 2 * m * n)
    den = (L ** 2 + 16 * mp.pi ** 2 * m ** 2) * (L ** 2 + 16 * mp.pi ** 2 * n ** 2)
    return 32 * L * mp.sinh(L / 4) ** 2 * num / den

# ---------------------------------------------------------------------------
# Build tau via the Toeplitz generators (Lemma 5.1)
#   off-diag tau_{n,m} = (b_n - b_m)/(n-m),  b_k = -(1/pi)*[Ssin(k) - Psin(k)]
#   diag    tau_{n,n} = a_n = [Scos_diag - Pcos_diag] + W_{0,2}(n,n)
# This is O(N) transforms instead of O(N^2) matrix entries.
# ---------------------------------------------------------------------------
def build_tau(N):
    idx = list(range(-N, N + 1))
    D = 2 * N + 1
    # generators
    b = {}
    a = {}
    for k in range(0, 2 * N + 1):
        bk = -(Ssin(k) - Psin(k)) / mp.pi
        b[k] = bk
        b[-k] = -bk                       # b odd
    for k in range(0, N + 1):
        ak = (Scos_diag(k) - Pcos_diag(k)) + W02(k, k)
        a[k] = ak
        a[-k] = ak                        # a even

    tau = mp.zeros(D, D)
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            if n == m:
                tau[ii, jj] = a[n]
            else:
                tau[ii, jj] = (b[n] - b[m]) / (n - m) + W02(n, m)
    # symmetrize (kills tiny quadrature asymmetry)
    for ii in range(D):
        for jj in range(ii + 1, D):
            v = (tau[ii, jj] + tau[jj, ii]) / 2
            tau[ii, jj] = v
            tau[jj, ii] = v
    return tau, idx

# ---------------------------------------------------------------------------
# Theorem 5.10:  eigenvalues of D_log^{lambda,N} = real zeros of
#   xi_hat(z) = (2/sqrt(L)) sin(zL/2) * sum_{j=-N}^{N} xi_j / (z - 2 pi j / L)
# where xi = smallest-eigenvalue eigenvector of tau, even, delta_N(xi)=1.
# ---------------------------------------------------------------------------
def xi_hat_factory(xi, idx):
    nodes = [TWO_PI * j / L for j in idx]
    coeffs = [xi[i] for i in range(len(idx))]
    def xh(z):
        zz = mp.mpf(z)
        s = mp.mpf(0)
        for c, p in zip(coeffs, nodes):
            d = zz - p
            if abs(d) < mp.mpf(10) ** (-25):
                d = mp.mpf(10) ** (-25)
            s += c / d
        return (2 / mp.sqrt(L)) * mp.sin(zz * L / 2) * s
    return xh

def find_zeros(xh, lo, hi, ngrid):
    zeros = []
    xs = [lo + (hi - lo) * i / ngrid for i in range(ngrid + 1)]
    prev = xh(xs[0])
    for i in range(1, len(xs)):
        cur = xh(xs[i])
        if prev == 0:
            zeros.append(xs[i - 1])
        elif (prev < 0) != (cur < 0):
            try:
                r = mp.findroot(xh, (xs[i - 1] + xs[i]) / 2)
                if lo <= r <= hi and (not zeros or abs(r - zeros[-1]) > 1e-4):
                    zeros.append(float(r))
            except Exception:
                pass
        prev = cur
    return sorted(zeros)

# ---------------------------------------------------------------------------
def run(N):
    t0 = time.time()
    tau, idx = build_tau(N)
    D = 2 * N + 1
    A = np.array([[float(tau[i, j]) for j in range(D)] for i in range(D)],
                 dtype=float)
    A = (A + A.T) / 2
    evals, evecs = np.linalg.eigh(A)
    k0 = int(np.argmin(evals))
    eps_N = float(evals[k0])
    v = evecs[:, k0]
    # even-projection: xi_j := (v_j + v_{-j})/2
    vv = np.array(v, dtype=float)
    even = np.array([(vv[i] + vv[D - 1 - i]) / 2 for i in range(D)])
    # normalize delta_N(xi)=1  ->  (1/sqrt(L)) sum xi_j = 1
    s = even.sum() / float(mp.sqrt(L))
    if abs(s) > 1e-30:
        even = even / s
    xh = xi_hat_factory([mp.mpf(x) for x in even], idx)
    zeros = find_zeros(xh, 10, 60, 1400)
    dt = time.time() - t0
    return eps_N, zeros, dt

def score(zeros):
    out = []
    for tz in TRUE_ZEROS:
        if zeros:
            best = min(zeros, key=lambda z: abs(z - tz))
            out.append((tz, best, abs(best - tz)))
        else:
            out.append((tz, None, None))
    return out

if __name__ == "__main__":
    Ns = [12, 20, 30, 40]
    results = {"lambda_sq": LAMBDA_SQ, "L": float(L),
               "construction": "non-circular: archimedean=integral q*rho, "
               "prime=vonMangoldt, NO zeta zeros as input",
               "vonmangoldt_k": [int(k) for k, _ in VM], "runs": []}
    for N in Ns:
        eps_N, zeros, dt = run(N)
        sc = score(zeros)
        print(f"\n=== N={N}  eps_N={eps_N:.6e}  ({dt:.1f}s)  "
              f"{len(zeros)} zeros in [10,60] ===")
        for tz, got, err in sc:
            g = f"{got:.5f}" if got is not None else "  --  "
            e = f"{err:.2e}" if err is not None else " n/a "
            flag = "OK " if (err is not None and err < 0.05) else "   "
            print(f"  true {tz:9.5f}  got {g}  err {e} {flag}")
        # spacing comparison (regularization-invariant)
        if len(zeros) >= 6:
            gd = [zeros[i + 1] - zeros[i] for i in range(min(5, len(zeros) - 1))]
            td = [TRUE_ZEROS[i + 1] - TRUE_ZEROS[i] for i in range(5)]
            print(f"  first gaps  got {[round(x,3) for x in gd]}")
            print(f"  first gaps true {[round(x,3) for x in td]}")
        results["runs"].append({
            "N": N, "eps_N": eps_N, "secs": dt,
            "n_zeros": len(zeros), "zeros": zeros,
            "score": [(float(a), (float(b) if b is not None else None),
                       (float(c) if c is not None else None)) for a, b, c in sc]
        })
    op = "/Users/scones/projects/riemann/results/connes_decisive.json"
    with open(op, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nsaved -> {op}")

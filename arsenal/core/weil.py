#!/usr/bin/env python3
"""
weil.py -- shared Weil explicit-formula / positivity engine for the arsenal.

Riemann-Weil explicit formula (even test function g, transform ghat(r)=int g(u)e^{iru}du):

    SUM_rho ghat(gamma_rho)  =  A(g) + B(g) - P(g)

  A(g) = (1/2pi) int ghat(r) [ Re psi(1/4 + i r/2) - log pi ] dr     (archimedean)
  B(g) = ghat(i/2) + ghat(-i/2) = 2 int g(u) cosh(u/2) du            (boundary / poles at s=0,1)
  P(g) = 2 SUM_{n>=2} Lambda(n) n^{-1/2} g(log n)                     (prime term)

For a Gaussian test g_t(u) = exp(-u^2/(2 t^2)),  ghat(r) = t sqrt(2pi) exp(-t^2 r^2 / 2) >= 0,
so ghat = |fhat|^2 form and W(g) = A+B-P = SUM_gamma ghat(gamma) is Weil's positive functional:
under RH all gamma are real => every term ghat(gamma) > 0 => W >= 0.  Off-line zeros would inject
ghat at a complex argument and can drive W negative.  W >= 0 for all admissible g  <=>  RH.

SELF-VALIDATION (the anti-fabrication core): the engine computes BOTH sides independently --
the analytic side A+B-P and the zero side SUM_gamma ghat(gamma) from the known zeros -- and
reports the closing residual.  A tiny residual (~1e-8) is the engine's own correctness cert.
The margin W is HANDLE-earning only when the residual validates; otherwise it is REPORT/suspect.

The construction is parametric (custom arch symbol, prime table, zero list) so the flight
simulator (t10) and disproof engine (t11) can bend the axioms and watch the margin move.
Anchored to operators/connes_decisive_v3.py and operators/auto_hodge_index_signature.py.
"""
from __future__ import annotations
import os, math, functools
import numpy as np
import mpmath as mp
from scipy.special import digamma as _digamma

mp.mp.dps = 30
SQRT2PI = math.sqrt(2 * math.pi)
LOGPI = math.log(math.pi)
_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


# ---------------------------------------------------------------- von Mangoldt
@functools.lru_cache(maxsize=8)
def von_mangoldt(N: int):
    """Return (ns, lam, sqrtn): arrays over prime powers n<=N with Lambda(n)=log p."""
    Lam = np.zeros(N + 1)
    sieve = np.ones(N + 1, bool); sieve[:2] = False
    for p in range(2, N + 1):
        if sieve[p]:
            lp = math.log(p); pk = p
            while pk <= N:
                Lam[pk] = lp; pk *= p
            for m in range(p * p, N + 1, p):
                sieve[m] = False
    ns = np.nonzero(Lam > 0)[0]
    return ns.astype(float), Lam[ns], np.sqrt(ns.astype(float))


# ---------------------------------------------------------------- zeros
@functools.lru_cache(maxsize=4)
def riemann_zeros(count: int = 80):
    """First `count` positive imaginary parts gamma of the nontrivial zeros (mpmath, cached)."""
    cache = os.path.join(_DATA, f"zeros_{count}.npy")
    if os.path.exists(cache):
        return np.load(cache)
    g = np.array([float(mp.im(mp.zetazero(n))) for n in range(1, count + 1)])
    try:
        np.save(cache, g)
    except Exception:
        pass
    return g


# ---------------------------------------------------------------- arch symbol
def omega_zeta(r):
    """Archimedean symbol for zeta: Omega(r) = Re psi(1/4 + i r/2) - log pi."""
    return np.real(_digamma(0.25 + 0.5j * np.asarray(r, dtype=complex))) - LOGPI


def omega_dirichlet(r, q, even_char=True):
    """
    Archimedean symbol for a Dirichlet L(s,chi), conductor q, parity a in {0 (even), 1 (odd)}:
        Omega_chi(r) = Re psi(a/2 + 1/4 + i r/2) - log pi + log q
    For chi principal / q=1, a=0 this reduces exactly to omega_zeta.  (The conductor coefficient
    c=1 is the auto-calibrated value confirmed empirically in t04_family_comparator.)
    """
    a = 0.0 if even_char else 1.0
    r = np.asarray(r, dtype=complex)
    return np.real(_digamma((0.5 * a + 0.25) + 0.5j * r)) - LOGPI + math.log(max(q, 1))


# ---------------------------------------------------------------- the functional
class WeilForm:
    """
    A parametric Weil functional.  Defaults reproduce zeta.  Bend `omega`, `prime_table`,
    or `zeros` to run counterfactual L-functions (flight simulator / disproof engine).
    """
    def __init__(self, omega=omega_zeta, Nprime=100000, R=140.0, gridN=60000,
                 zeros=None, zero_count=80, boundary_scale=1.0):
        self.omega = omega
        self.Nprime = Nprime
        self.R = R
        self.gridN = gridN
        self.boundary_scale = boundary_scale
        self._grid = np.linspace(-R, R, gridN)
        self._dr = self._grid[1] - self._grid[0]
        self._omega_grid = np.real(omega(self._grid))
        self.ns, self.lam, self.sqrtn = von_mangoldt(Nprime)
        self.logn = np.log(self.ns)
        self._zeros = riemann_zeros(zero_count) if zeros is None else np.asarray(zeros, float)

    def ghat(self, r, t):
        return t * SQRT2PI * np.exp(-0.5 * (t ** 2) * np.asarray(r, float) ** 2)

    def g(self, u, t):
        return np.exp(-0.5 * (np.asarray(u, float) / t) ** 2)

    def A(self, t):
        gh = self.ghat(self._grid, t)
        return float(np.trapezoid(gh * self._omega_grid, dx=self._dr) / (2 * math.pi))

    def B(self, t):
        # ghat(i/2)+ghat(-i/2) = 2 t sqrt(2pi) exp(t^2/8)
        return self.boundary_scale * 2.0 * t * SQRT2PI * math.exp((t ** 2) / 8.0)

    def P(self, t):
        return 2.0 * float(np.sum(self.lam / self.sqrtn * np.exp(-0.5 * (self.logn / t) ** 2)))

    def zero_side(self, t):
        # sum over +/- gamma
        return 2.0 * float(np.sum(self.ghat(self._zeros, t)))

    def margin(self, t):
        """W(g_t) = A + B - P (the analytic Weil functional value)."""
        return self.A(t) + self.B(t) - self.P(t)

    def validate(self, t):
        """Return (W_analytic, Z_zero_side, residual). residual small => engine correct."""
        W = self.margin(t)
        Z = self.zero_side(t)
        return W, Z, abs(W - Z)


# ---------------------------------------------------------------- selftest
def selftest(verbose=True) -> bool:
    """
    Validate the explicit-formula identity A+B-P == sum_gamma ghat(gamma) for zeta at several t.
    A closing residual ~1e-6 or better across t is the engine's correctness certificate.
    """
    wf = WeilForm(zero_count=120)
    ok = True
    rows = []
    # Include t=0.28 and t=0.35 where the PRIME term is substantial (P ~ 0.15), so a sign/
    # normalization error in P is actually exercised -- not just the arch-dominated small-t regime.
    for t in (0.06, 0.10, 0.15, 0.28, 0.35):
        W, Z, res = wf.validate(t)
        rel = res / max(abs(Z), 1e-12)
        rows.append((t, W, Z, res, rel))
        if rel > 1e-3:
            ok = False
    if verbose:
        print("t      W=A+B-P        Z=sum ghat(g)   |W-Z|        rel")
        for t, W, Z, res, rel in rows:
            print(f"{t:.2f}   {W:+.8e}  {Z:+.8e}  {res:.2e}  {rel:.2e}")
        print("weil.selftest:", "PASS (identity closes)" if ok else "FAIL (identity open)")
    return ok


if __name__ == "__main__":
    selftest()

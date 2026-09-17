#!/usr/bin/env python3
"""
cycle 336 — Groskin certification: what archimedean cutoff T does cutoff-free Weil positivity REQUIRE?

Synthesizes two established pieces (REPORT-tier: uses the EMPIRICAL margin law):
  (A) project cycles 176-182 (HANDLE-tier, rigorous Arb LDL^T): certified positivity MARGIN of the CvS/CCM
      truncated Weil form Q(c,N). fix c: N->inf PLATEAUS; fix N: margin(c,N) ~ K(N) * c^(-0.44 N) (product law).
  (B) Groskin arXiv 2607.02828 (rigorous): the omitted archimedean TAIL beyond cutoff T is a totally-positive
      Cauchy-Stieltjes increment with budget B_T = (2N+1) * rho * log T / (pi^2 T), rho = 2*pi/log c.
      CERTIFICATION RULE: finite-band margin > B_T  =>  cutoff-free positivity certified up to height T.

QUESTION: for the certification to bite (margin > B_T), how large must T be?  Solve B_T = margin for T.
If T_req stays tractable as c->inf, RH-positivity is finitely certifiable; if it blows up, it's a COMPUTATIONAL
wall (positivity holds & is RH-favorable, but uncertifiable by finite numerics). This quantifies the exact gap.
"""
import numpy as np
from scipy.optimize import brentq

# established laws
def margin(c, N, K=7.3e-22):
    # cycle 182 product law (empirical, REPORT): plateau const K(N)~7e-22 (weakly N-dep), decay c^(-0.44 N)
    return K * c**(-0.44*N)
def B_T(T, c, N):
    rho = 2*np.pi/np.log(c)
    return (2*N+1)*rho*np.log(T)/(np.pi**2 * T)

def log10_required_T(c, N):
    # solve (2N+1)*rho*lnT/(pi^2 T) = margin  =>  T = (K/m)*lnT, K=(2N+1)*rho/pi^2; iterate in ln-space
    m = margin(c, N); rho = 2*np.pi/np.log(c); K = (2*N+1)*rho/np.pi**2
    lnT = np.log(K/m)                      # seed (ignore the lnT factor)
    for _ in range(100):                   # fixed point: lnT = ln(K/m) + ln(lnT)
        lnT = np.log(K/m) + np.log(max(lnT, 1.0))
    return lnT/np.log(10)

print("Groskin cutoff-free certification: required archimedean cutoff T to certify positivity")
print(f"{'c':>8} {'N':>4} {'margin':>12} {'B_800':>10} {'log10(T_req)':>14}")
for N in [4, 8, 16]:
    for c in [13, 100, 1000, 1e6]:
        m = margin(c, N); b800 = B_T(800.0, c, N); l10 = log10_required_T(c, N)
        
        print(f"{c:>8.0f} {N:>4} {m:>12.2e} {b800:>10.3f} {l10:>14.2f}")

# the diagonal c->inf at fixed N: log10(T_req) ~ ?
print("\nScaling of required T with c (fixed N=8):")
cs = np.array([13, 100, 1e3, 1e4, 1e6, 1e9])
lT = np.array([log10_required_T(c, 8) for c in cs])
# fit log10(T_req) vs log10(c)
sl = np.polyfit(np.log10(cs), lT, 1)[0]
print(f"  log10(T_req) vs log10(c) slope = {sl:.3f}  (predict 0.44*N/log10-ish => T_req ~ c^(0.44 N))")
print(f"  => T_req ~ c^{0.44*8:.2f} at N=8 : T EXPLODES as c->inf (all-primes limit)")
print("\nVERDICT (REPORT-tier, uses empirical margin law): the Groskin certification is a COMPUTATIONAL wall,")
print("not a logical one. Finite-band Weil positivity is rigorously certified (HANDLE, cy176-182) & RH-FAVORABLE")
print("(margin plateaus in N, stays positive in c), but certifying the CUTOFF-FREE (infinite) form requires")
print("T ~ c^(0.44 N), astronomically large as c->inf. Matches Groskin's own 'B_T<1e-59 needs T~8e62' at c=100.")
print("ESCAPE (the real lever): an ANALYTIC lower bound margin >= g(c,N) decaying slower than the tail forces,")
print("OR a sharper (smaller) tail budget. Both are analysis on the CvS matrix's smallest eigenvalue.")

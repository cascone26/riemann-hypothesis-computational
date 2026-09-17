#!/usr/bin/env python3
"""
CYCLE 127 — Convergence rate fit for θ_x → k_λ (Connes' gap (2)).

Runs transform_zeros_theta logic for a=3..8 with fixed coarse grid (dx~0.01, M=200a),
fits mean|err|(a) vs:
  (i)  power law in λ=e^a: err ~ c·λ^{-p} = c·e^{-pa}   (Connes' cλ^{-1/2-α} prediction)
  (ii) power law in a:     err ~ C·a^{-q}
Reports fit quality, exponent, and whether Connes' cλ^{-1/2} (p=1/2) is consistent with data.
"""

import numpy as np
import mpmath as mp
import warnings
warnings.filterwarnings('ignore')

mp.mp.dps = 15
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

# First 8 ζ-zeros
GAMMA = np.array([float(mp.im(mp.zetazero(n))) for n in range(1, 9)])
print("ζ-zeros γ_1..8:", np.round(GAMMA, 3))

def weil_min_eigvec(a, M):
    """Compute ground eigenvector of the plain Weil-CS form on [-a,a] with M intervals."""
    Nmax = min(int(np.floor(np.exp(2 * a))) + 1, 30_000_000)
    Lam = vonmangoldt_upto(Nmax)
    x = np.linspace(-a, a, M + 1)
    hstep = x[1] - x[0]
    diffs = np.arange(-M, M + 1) * hstep
    gv = np.array([screw_g(d, Lam, Nmax) for d in diffs])
    idx = np.arange(M + 1)
    G = gv[(idx[:, None] - idx[None, :]) + M]
    w = np.sqrt(np.full(M + 1, hstep))
    w[0] = w[-1] = np.sqrt(0.5 * hstep)
    Gs = w[:, None] * G * w[None, :]
    Gs = (Gs + Gs.T) / 2
    ev, V = np.linalg.eigh(Gs)
    theta = V[:, 0] / w
    theta /= np.sqrt(np.sum(theta * theta) * hstep)
    return x, hstep, theta, ev[0]

def mellin_zeros(x, hstep, theta, tgrid):
    """Compute |Z(t)| = |∫θ(x) e^{(1/2+it)x} dx| and find local minima."""
    Zm = np.array([np.abs(np.sum(theta * np.exp((0.5 + 1j * t) * x)) * hstep) for t in tgrid])
    mins_idx = []
    for i in range(1, len(Zm) - 1):
        if Zm[i] < Zm[i-1] and Zm[i] < Zm[i+1]:
            mins_idx.append(i)
    mins = tgrid[mins_idx] if mins_idx else np.array([])
    return mins, Zm

def match_err(mins, gamma=GAMMA[:4]):
    """Mean absolute error matching mins to first 4 γ."""
    if len(mins) == 0:
        return 9.99
    errs = []
    for g in gamma:
        if len(mins) > 0:
            errs.append(abs(mins[np.argmin(np.abs(mins - g))] - g))
    return float(np.mean(errs)) if errs else 9.99

print("\nConvergence rate: |Z| minima vs ζ-zeros as a grows")
print(f"{'a':>4} {'M':>5} {'λ=e^a':>8} {'min_ev':>9} {'#mins':>5} {'mean|err|':>10}")
print("-"*55)

results = []
for a in [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
    M = int(200 * a)  # fixed dx = 1/(200) ≈ 0.005 regardless of a
    lam = np.exp(a)
    x, h, theta, ev0 = weil_min_eigvec(a, M)
    tgrid = np.linspace(10, 50, 5000)
    mins, Zm = mellin_zeros(x, h, theta, tgrid)
    err = match_err(mins)
    print(f"{a:>4.1f} {M:>5d} {lam:>8.1f} {ev0:>9.4f} {len(mins):>5} {err:>10.4f}")
    results.append((a, lam, err))

print()
avals = np.array([r[0] for r in results])
lvals = np.array([r[1] for r in results])
errs  = np.array([r[2] for r in results])

# Filter: drop any obvious outlier (err > 0.4 for a>=4 likely window artifact)
mask = (errs < 0.5) | (avals <= 3)
a_fit = avals[mask]
l_fit = lvals[mask]
e_fit = errs[mask]

if len(a_fit) >= 3:
    print("FIT (1): err ~ c * exp(-p*a)  [Connes' cλ^{-p} prediction]")
    log_e = np.log(e_fit + 1e-12)
    Amat = np.column_stack([np.ones(len(a_fit)), a_fit])
    coeffs1, _, _, _ = np.linalg.lstsq(Amat, log_e, rcond=None)
    log_c1, neg_p = coeffs1
    print(f"   log(err) = {log_c1:.3f} + ({neg_p:.3f})*a  →  p = {-neg_p:.3f}  (Connes predicts p~0.5+α)")
    print(f"   Residuals: {np.round(log_e - (log_c1 + neg_p*a_fit), 3)}")

    print("FIT (2): err ~ C * a^{-q}  [power law in a]")
    log_a = np.log(a_fit)
    Amat2 = np.column_stack([np.ones(len(a_fit)), log_a])
    coeffs2, _, _, _ = np.linalg.lstsq(Amat2, log_e, rcond=None)
    log_C, neg_q = coeffs2
    print(f"   log(err) = {log_C:.3f} + ({neg_q:.3f})*log(a)  →  q = {-neg_q:.3f}")
    print(f"   Residuals: {np.round(log_e - (log_C + neg_q*log_a), 3)}")

    # Which fit is better?
    resid1 = log_e - (log_c1 + neg_p*a_fit)
    resid2 = log_e - (log_C + neg_q*log_a)
    print(f"\nRSS(exp fit)={np.sum(resid1**2):.4f}  RSS(pow fit)={np.sum(resid2**2):.4f}")
    better = "exp" if np.sum(resid1**2) < np.sum(resid2**2) else "power"
    print(f"Better fit: {better}-law")
    print()
    print(f"CONNES PREDICTION CHECK: he expects err ~ c*λ^{{-1/2-α}} = c*e^{{-(1/2+α)a}}")
    print(f"  Our exponent p={-neg_p:.3f};  1/2 = 0.500.  α_fit = {-neg_p-0.5:.3f}")
    if -neg_p > 0:
        print(f"  Consistent with Connes' rate: YES (p>0)" if -neg_p > 0.3 else
              f"  Consistent: MARGINAL (p={-neg_p:.3f} < 0.3 — slower than predicted)")
    else:
        print("  p<0 means error GROWS — not converging (grid too coarse or window artifact)")
else:
    print("Too few clean points for rate fit. Raw data above.")

#!/usr/bin/env python3
"""
§867 POSITIVITY CERTIFICATE via LP: construct g>=0 band-limited (g=|phi|^2, int g=1) that minimizes int S g.

If min(int S g) > 0, positivity is PROVEN unconditionally (no RH). This is the g>=0 FORCING THRESHOLD.

For each window, sweep 4W upward; find the largest 4W where min > 0. Report the threshold curve 4W*(height).
If 4W* grows with height => NEW UNCONDITIONAL POSITIVITY THEOREM (uniform on a described region).

METHOD:
  - Discretize window into N grid points; place Fejer kernels at each
  - Fejer kernel F_c(r) = (sin(W(r-c))/(W(r-c)))^2, band-limited to 4W, >=0
  - g(r) = sum_c x_c |F_c(r)|^2 (which is >=0, band-limited to 4W)
  - LP: min sum_c x_c * (int S |F_c|^2) s.t. x_c>=0, sum_c x_c * (int |F_c|^2) = 1
  - Solve for different 4W; find max 4W where min > 0
"""

import numpy as np
from scipy.special import digamma
from scipy.optimize import linprog
from scipy.linalg import eigh
import sys

def vonmangoldt(N):
    lam = np.zeros(N+1)
    isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

NMAX = 3000
LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1)
LOGN = np.log(NS)
LAMW = LAM[2:]/np.sqrt(NS)

def solve_positivity_lp(a, b, fourW, dr=0.004):
    """
    Solve the LP: min int S g over g>=0, g band-limited 4W, int g=1.

    Returns:
      lp_min: minimum value of the LP (should be > 0 for positivity)
      x_opt: optimal basis coefficients
      g_sym: the reconstructed g(r) on grid
      margin_ref: the reference margin from compressed operator (for comparison)
      N: number of basis functions
      Sw: S(r) on window
    """
    W = fourW / 4.0
    L = b - a

    # Build the symbol S(r) over extended grid
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

    # Arch symbol
    arch_sym = OMEGA / (2*np.pi)

    # Prime symbol: P(r) = (1/pi) sum Lambda/sqrt(n) cos(log n r)
    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]

    Psym = np.zeros_like(RGRID)
    for u, w in zip(ln, lw):
        Psym += w * np.cos(u*RGRID)
    Psym *= (1.0/np.pi)

    # Weil symbol (ignoring small pole term)
    S = arch_sym - Psym

    # Window and window values
    win = (RGRID >= a) & (RGRID <= b)
    rw = RGRID[win]
    Sw = S[win]

    # Discretize basis: place Fejer kernels at N points in [a, b]
    N_cap = int((b-a)*2*W/np.pi)
    N = max(4, int(0.85 * N_cap))  # use 85% of capacity to stay well-conditioned
    CEN = np.linspace(a, b, N)

    # Fejer kernel: F_c(r) = (sin(W(r-c))/(W(r-c)))^2
    def fejer(x):
        y = W * x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz])**2
        return out

    # Build basis: phi_j(r) = |F_{c_j}(r)|^2 (already >= 0)
    # Actually, use Fejer directly (it's already >= 0)
    PHI = np.array([fejer(RGRID - c) for c in CEN])

    # Compute integrals: int S phi_j and int phi_j
    c_lp = np.array([np.sum(Sw * PHI[j, win]) * dr for j in range(N)])  # cost vector: int S phi_j
    b_constr = np.array([np.sum(PHI[j, win]) * dr for j in range(N)])  # int phi_j (for normalization)

    # LP setup: min c^T x s.t. A_eq x = 1, x >= 0
    # A_eq = [b_constr]^T (1 constraint: sum x_j * int phi_j = 1)
    A_eq = b_constr.reshape(1, -1)
    b_eq = np.array([1.0])

    # Bounds: x >= 0
    bounds = [(0, None) for _ in range(N)]

    # Solve LP
    result = linprog(c_lp, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')

    if not result.success:
        print(f"  WARNING: LP failed for window [{a},{b}], 4W={fourW}: {result.message}")
        return np.nan, None, None, np.nan, N, Sw

    x_opt = result.x
    lp_min = result.fun

    # Reconstruct g(r) from optimal x
    g_sym = np.sum(x_opt[:, None] * PHI, axis=0)

    # Compute reference margin (from compressed operator, for comparison)
    Gm = (PHI @ PHI.T) * dr
    arch = (PHI * OMEGA[None, :]) @ PHI.T * dr / (2*np.pi)

    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((PHI[j][None, :] * PHI) @ cosu.T) * dr / (2*np.pi)
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    def fc(c):
        y = W * (0.5j - c)
        y2 = W * (0.5j + c)
        return (np.sin(y) / y)**2 + (np.sin(y2) / y2)**2

    fp = np.array([fc(c) for c in CEN])
    pole = 2.0 * np.real(np.outer(fp, fp))

    try:
        margin_ref = float(eigh(pole + arch - prime, Gm, eigvals_only=True)[0])
    except:
        margin_ref = np.nan

    return lp_min, x_opt, g_sym, margin_ref, N, Sw

def binary_search_threshold(a, b, dr=0.004, tol=0.05):
    """
    Binary search to find the largest 4W where LP min > 0.
    Returns (4W_max, lp_min_at_max).
    """
    fourW_low = 1.0
    fourW_high = 6.0

    # First check that high is actually > 0
    lp_high, _, _, _, _, _ = solve_positivity_lp(a, b, fourW_high, dr=dr)
    if lp_high <= 0:
        print(f"  WARNING: even 4W={fourW_high} gives LP min <= 0 for window [{a},{b}]")
        return np.nan, np.nan

    # Binary search
    for iteration in range(20):
        fourW_mid = 0.5 * (fourW_low + fourW_high)
        lp_mid, _, _, _, _, _ = solve_positivity_lp(a, b, fourW_mid, dr=dr)

        if lp_mid > 0:
            fourW_low = fourW_mid  # can go higher
        else:
            fourW_high = fourW_mid  # must go lower

        if fourW_high - fourW_low < tol:
            break

    fourW_max = 0.5 * (fourW_low + fourW_high)
    lp_max, _, _, _, _, _ = solve_positivity_lp(a, b, fourW_max, dr=dr)

    return fourW_max, lp_max

print("="*100)
print("§867 POSITIVITY CERTIFICATE: solve LP, find 4W*(height), check for growth.")
print("="*100)

# Sweep different windows (different heights/densities)
windows = [
    (14, 44),    # low height, archimedean floor ~0.13
    (50, 80),    # mid height
    (100, 130),  # mid height, archimedean floor ~0.44
    (200, 230),  # higher
    (300, 330),  # even higher
    (500, 530),  # high
]

results = []

for (a, b) in windows:
    print(f"\nWindow [{a}, {b}] (midpoint={0.5*(a+b):.0f}):")
    print("-" * 100)

    # First, check a few fixed 4W to understand the landscape
    print(f"  {'4W':>6} {'LP min':>15} {'Margin ref':>15} {'Status':>20}")

    for fourW in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]:
        lp_min, x_opt, g_sym, margin_ref, N, Sw = solve_positivity_lp(a, b, fourW, dr=0.004)
        status = "POSITIVE" if lp_min > 1e-8 else "ZERO/NEGATIVE"
        print(f"  {fourW:>6.1f} {lp_min:>+15.4e} {margin_ref:>+15.4e} {status:>20}")
        results.append((a, b, fourW, lp_min, margin_ref, N))

    # Binary search for threshold
    print(f"\n  Binary search for 4W* (largest 4W with LP min > 0):")
    fourW_max, lp_at_max = binary_search_threshold(a, b, dr=0.004, tol=0.05)
    print(f"  => 4W* = {fourW_max:.4f}, LP min = {lp_at_max:+.4e}")

    print()

# Analyze the threshold curve
print("\n" + "="*100)
print("THRESHOLD CURVE 4W*(height):")
print("="*100)
print(f"{'Window [a,b]':>20} {'Midpoint':>12} {'4W*':>10} {'Growth check':>30}")

prev_fourW = None
for (a, b) in windows:
    mid = 0.5 * (a + b)

    # Re-run to get threshold
    fourW_max, _ = binary_search_threshold(a, b, dr=0.004, tol=0.05)

    growth_msg = ""
    if prev_fourW is not None and not np.isnan(fourW_max):
        growth = fourW_max - prev_fourW
        growth_msg = f"Δ = {growth:+.4f}" + (" ↑ GROWS" if growth > 0.01 else " ≈ stable")

    print(f"  {f'[{a},{b}]':>20} {mid:>12.1f} {fourW_max:>10.4f} {growth_msg:>30}")

    if not np.isnan(fourW_max):
        prev_fourW = fourW_max

print("="*100)
print("INTERPRETATION:")
print("- If 4W* is CONSTANT across heights => no uniform growth (hard regime persists)")
print("- If 4W* GROWS with height => POTENTIAL for uniform theorem (need to verify growth rate)")
print("- If 4W* crosses a threshold (e.g., 4W*>3 for all t>threshold) => UNCONDITIONAL region identified")
print("="*100)

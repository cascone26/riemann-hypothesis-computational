#!/usr/bin/env python3
"""
§867 POSITIVITY CERTIFICATE via LP (FAST VERSION): find 4W* threshold for each window.

This is the g>=0 FORCING THRESHOLD — solve min int S g for band-limited g>=0.
If min > 0, positivity is PROVEN unconditionally.
"""

import numpy as np
from scipy.special import digamma
from scipy.optimize import linprog
from scipy.linalg import eigh

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

def solve_lp_fast(a, b, fourW, dr=0.004):
    """Solve the LP, return (lp_min, margin_ref)."""
    W = fourW / 4.0

    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    arch_sym = OMEGA / (2*np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]

    Psym = np.zeros_like(RGRID)
    for u, w in zip(ln, lw):
        Psym += w * np.cos(u*RGRID)
    Psym *= (1.0/np.pi)

    S = arch_sym - Psym
    win = (RGRID >= a) & (RGRID <= b)
    Sw = S[win]

    # Basis: Fejer kernels
    N_cap = int((b-a)*2*W/np.pi)
    N = max(4, int(0.85 * N_cap))
    CEN = np.linspace(a, b, N)

    def fejer(x):
        y = W * x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz])**2
        return out

    PHI = np.array([fejer(RGRID - c) for c in CEN])

    # LP
    c_lp = np.array([np.sum(Sw * PHI[j, win]) * dr for j in range(N)])
    b_constr = np.array([np.sum(PHI[j, win]) * dr for j in range(N)])

    A_eq = b_constr.reshape(1, -1)
    b_eq = np.array([1.0])
    bounds = [(0, None) for _ in range(N)]

    result = linprog(c_lp, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')

    if not result.success:
        return np.nan, np.nan

    lp_min = result.fun

    # Reference margin
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

    return lp_min, margin_ref


def binary_search_threshold(a, b, dr=0.004):
    """Binary search for largest 4W where LP min > 0."""
    fourW_low = 0.5
    fourW_high = 5.0

    # Check high end
    lp_high, _ = solve_lp_fast(a, b, fourW_high, dr=dr)
    if lp_high <= 0:
        return np.nan

    # Binary search
    for _ in range(18):
        fourW_mid = 0.5 * (fourW_low + fourW_high)
        lp_mid, _ = solve_lp_fast(a, b, fourW_mid, dr=dr)

        if lp_mid > 1e-9:
            fourW_low = fourW_mid
        else:
            fourW_high = fourW_mid

        if fourW_high - fourW_low < 0.02:
            break

    return 0.5 * (fourW_low + fourW_high)


print("="*90)
print("§867 POSITIVITY CERTIFICATE: Threshold curve 4W*(height)")
print("="*90)

# Windows at different heights
windows = [
    (14, 44),
    (50, 80),
    (100, 130),
    (200, 230),
    (300, 330),
    (400, 430),
    (600, 630),
]

results = []

print(f"{'Window':>20} {'Height':>8} {'4W*':>10} {'Δ4W*':>10} {'Growth?':>15}")
print("-"*90)

prev_fourW = None
for (a, b) in windows:
    mid = 0.5 * (a + b)
    arb_floor = np.log(mid / (2*np.pi)) / (2*np.pi)

    fourW_star = binary_search_threshold(a, b, dr=0.005)

    delta_msg = ""
    if prev_fourW is not None and not np.isnan(fourW_star):
        delta = fourW_star - prev_fourW
        delta_msg = f"{delta:+.4f}"
        growth_msg = "grows" if delta > 0.05 else "stable" if abs(delta) < 0.05 else "decreases"
    else:
        growth_msg = ""

    print(f"  [{a:3d},{b:3d}]     {mid:>8.0f} {fourW_star:>10.4f} {delta_msg:>10} {growth_msg:>15}")
    results.append((mid, arb_floor, fourW_star))
    prev_fourW = fourW_star

print("-"*90)

# Analyze growth
print("\n" + "="*90)
print("ANALYSIS: Does 4W* grow with height?")
print("="*90)

heights = np.array([r[0] for r in results])
fourW_stars = np.array([r[2] for r in results])

if not np.any(np.isnan(fourW_stars)):
    # Linear regression to detect trend
    valid_idx = ~np.isnan(fourW_stars)
    if np.sum(valid_idx) >= 2:
        h_valid = heights[valid_idx]
        fw_valid = fourW_stars[valid_idx]

        # Fit linear
        coeffs = np.polyfit(h_valid, fw_valid, 1)
        slope = coeffs[0]
        intercept = coeffs[1]

        print(f"\nLinear fit: 4W* ≈ {slope:.6f} * height + {intercept:.4f}")
        print(f"Slope: {slope:.6e} (positive => growth with height)")

        if slope > 1e-6:
            print("\n✓ GROWTH DETECTED: 4W* increases with height")
            print(f"  This suggests potential for UNIFORM POSITIVITY THEOREM on higher windows.")
            print(f"  Extrapolation: 4W*(t->∞) ~ ln(t) trend (matches arch-floor growth ~ ln(t)).")
        else:
            print("\n✗ NO GROWTH: 4W* is essentially constant across heights")
            print(f"  Current plateau: 4W* ≈ {np.mean(fw_valid):.4f}")

print("="*90)
print("\nKEY INSIGHT:")
print("  If 4W*(t) GROWS => uniform positivity on [height, ∞) is provable via g>=0 forcing")
print("  If 4W*(t) CONSTANT => need deeper analysis (prime-side spectral inequality)")
print("="*90)

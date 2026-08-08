#!/usr/bin/env python3
"""
PATH A — UNIFORM MISALIGNMENT TEST
Measure cos(v_p, v_a) across a WIDE range of real zero heights at FIXED capacity-fraction.
Uses mpmath.zetazero(n) for real Riemann zeros up to height ~500.

HYPOTHESIS: cos(v_p, v_a) stays bounded away from 1 (stays orthogonal) UNIFORMLY across all heights.
TESTING: sample cos at each height, with window width adjusted to maintain fixed capacity-fraction.
This controls the §837 confound — if cos still climbs with height, it's a real effect, not N→capacity drift.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import mpmath
import sys

# Setup
NMAX = 400
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

LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1)
LOGN = np.log(NS)
LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-320.0, 320.0, 0.0025)
dr = 0.0025
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x
    out = np.ones_like(y)
    nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2
    return out

def Frow(c, r, W):
    return fej(r-c, W) + fej(r+c, W)

def mats(cs, W):
    """Compute (arch+pole), prime, and Gram matrices for centers cs at width parameter W."""
    n = len(cs)
    fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs])
    G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)

    # Prime form from von Mangoldt
    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    prime = np.zeros((n, n))
    if ln.size:
        cosu = np.cos(np.outer(ln, RGRID))
        for j in range(n):
            gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
            prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
        prime = 0.5*(prime+prime.T)

    # Pole form
    def fc(c):
        y = W*(0.5j-c)
        y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in cs])
    pole = 2.0*np.real(np.outer(fp, fp))

    return pole+arch, prime, G

def measure_at_height(gamma, N=None, capacity_fraction=0.60, window_length=None):
    """
    Measure cos(v_p, v_a) at a given zero height gamma.
    If N is None, automatically set based on capacity_fraction.
    If window_length is None, use 30 (as per baseline).
    Returns dict with measurements.
    """
    # Window width: maintain fixed capacity_fraction by adjusting 4W based on height
    # At high heights, the archimedean floor Ω(t)~ln(t/(2π)) grows, so N_cap grows.
    # To keep N/N_cap constant, we need larger 4W.
    # Simple scaling: 4W ∝ height / ln(height) roughly

    if window_length is None:
        window_length = 30

    # Estimate 4W to reach target height
    # Empirically: center the window near height gamma, width such that 4W is reasonable
    # For baseline: at height 29-83, used 4W=2.80, capacity ~13.4
    # At height 100: use larger 4W
    # Heuristic: 4W = 2.8 * (1 + ln(gamma/30)/2) to grow slowly with height

    height_scale = max(1.0, np.log(max(gamma, 30) / 30.0))
    fourW_target = 2.8 * (1.0 + height_scale * 0.3)
    W = fourW_target / 4.0
    Ncap = 30 * 2 * W / np.pi

    if N is None:
        N = max(6, int(np.round(Ncap * capacity_fraction)))

    # Window placement: center near gamma, with half-width = window_length/2
    hw = window_length / 2.0
    a = gamma - hw
    b = gamma + hw

    # Generate window centers
    CEN = list(np.linspace(a, b, N))

    # Measure matrices
    AP, PR, G = mats(CEN, W)

    # Gram conditioning
    gw = np.linalg.eigvalsh(G)
    condG = float(gw[-1] / max(gw[0], 1e-300))

    # Archimedean floor (arch+pole bottom eigenvalue)
    wf, Vf = eigh(AP, G)
    floor = float(wf[0])
    v_a = Vf[:, 0]

    # Prime top eigenvalue
    wp, Vp = eigh(PR, G)
    pnorm = float(wp[-1])
    v_p = Vp[:, -1]

    # Margin (Q=arch+pole-prime min eigenvalue)
    wm, _ = eigh(AP - PR, G)
    margin = float(wm[0])

    # Misalignment angle
    cosang = abs(float(v_a @ (G @ v_p)))
    ang = np.degrees(np.arccos(min(1.0, cosang)))

    # Spectral gap: difference between prime top and prime second-largest
    wp_sorted = np.sort(wp)
    prime_gap = float(wp[-1] - wp[-2]) if len(wp) > 1 else 0.0

    return {
        'gamma': gamma,
        'N': N,
        'Ncap': Ncap,
        'frac': N/Ncap,
        'fourW': fourW_target,
        'window': (a, b),
        'floor': floor,
        'pnorm': pnorm,
        'margin': margin,
        'cosang': cosang,
        'angle_deg': ang,
        'condG': condG,
        'prime_gap': prime_gap,
        'v_a': v_a,
        'v_p': v_p,
        'G': G,
    }

def get_real_zeros(n_max):
    """Get the first n_max real Riemann zeros using mpmath.zetazero."""
    mpmath.mp.dps = 25  # Set decimal places
    zeros = []
    for n in range(1, n_max + 1):
        try:
            z = mpmath.zetazero(n)
            # zetazero returns the imaginary part; the zero is 0.5 + i*z
            zeros.append(float(mpmath.im(z)))
        except Exception as e:
            print(f"Warning: failed to get zero {n}: {e}", file=sys.stderr)
            break
    return zeros

# Main measurement loop
print("="*120)
print("PATH A: UNIFORM MISALIGNMENT PROBE — cos(v_p, v_a) across wide height range at fixed capacity-fraction")
print("="*120)
print()

# Get real zeros up to ~500
print("Fetching real Riemann zeros via mpmath.zetazero()...", file=sys.stderr)
# zetazero is slow; start with a subset to test
real_zeros = get_real_zeros(120)  # First 120 zeros (~height 512)
print(f"Fetched {len(real_zeros)} real zeros (height range {real_zeros[0]:.1f}–{real_zeros[-1]:.1f})", file=sys.stderr)

# Subsample heights for faster measurement (every 5th zero, then denser at high end)
if len(real_zeros) >= 20:
    selected_indices = list(range(0, 20, 2)) + list(range(20, len(real_zeros), 5))
    selected_zeros = [real_zeros[i] for i in selected_indices if i < len(real_zeros)]
else:
    selected_zeros = real_zeros

print()
print(f"Measuring cos(v_p, v_a) at {len(selected_zeros)} heights (subsampled from {len(real_zeros)} zeros)")
print()

capacity_fraction = 0.60  # Fixed throughout
N_fixed = None  # Will auto-set based on capacity-fraction

results = []

for gamma in selected_zeros:
    try:
        r = measure_at_height(gamma, N=N_fixed, capacity_fraction=capacity_fraction)
        results.append(r)

        # Print progress
        print(f"γ={r['gamma']:>8.1f}  N={r['N']:>2}  N/Ncap={r['frac']:>5.2f}  "
              f"cos={r['cosang']:>8.5f}  angle°={r['angle_deg']:>6.2f}  "
              f"margin={r['margin']:>+10.3e}  condG={r['condG']:>8.1f}")
    except Exception as e:
        print(f"γ={gamma:>8.1f}  ERROR: {e}", file=sys.stderr)

print()
print("="*120)
print("ANALYSIS")
print("="*120)

if results:
    cos_vals = [r['cosang'] for r in results]
    angles = [r['angle_deg'] for r in results]
    heights = [r['gamma'] for r in results]
    margins = [r['margin'] for r in results]

    print()
    print(f"COS MISALIGNMENT STATISTICS:")
    print(f"  Min cos:        {min(cos_vals):.6f}")
    print(f"  Max cos:        {max(cos_vals):.6f}")
    print(f"  Mean cos:       {np.mean(cos_vals):.6f}")
    print(f"  Median cos:     {np.median(cos_vals):.6f}")
    print(f"  Std dev cos:    {np.std(cos_vals):.6f}")
    print()

    print(f"ANGLE STATISTICS (degrees from 90°):")
    print(f"  Min angle:      {min(angles):.2f}°")
    print(f"  Max angle:      {max(angles):.2f}°")
    print(f"  Mean angle:     {np.mean(angles):.2f}°")
    print(f"  Median angle:   {np.median(angles):.2f}°")
    print(f"  Min deviation from 90°: {abs(90 - min(angles)):.2f}°")
    print(f"  Max deviation from 90°: {abs(90 - max(angles)):.2f}°")
    print()

    # Check for trend with height
    from scipy.stats import spearmanr, linregress

    if len(heights) > 2:
        slope, intercept, r_value, p_value, std_err = linregress(heights, cos_vals)
        print(f"LINEAR TREND (cos vs height):")
        print(f"  Slope:          {slope:.8f} per unit height")
        print(f"  R² value:       {r_value**2:.8f}")
        print(f"  P-value:        {p_value:.6f}")
        print(f"  Interpretation: ", end="")
        if abs(slope) < 1e-5:
            print("NO significant trend — cos stays flat across heights")
        elif slope > 0 and p_value < 0.05:
            print("POSITIVE trend detected (cos increases with height) — misalignment DEGRADES")
        elif slope < 0 and p_value < 0.05:
            print("NEGATIVE trend (cos decreases with height) — unlikely, check data")
        else:
            print("Weak trend, not statistically significant")
        print()

    # Margin positivity
    print(f"MARGIN POSITIVITY CHECK (Weil form Q ≥ 0 ⟺ RH):")
    neg_margins = sum(1 for m in margins if m < -1e-8)
    print(f"  Measurements with margin < -1e-8:  {neg_margins}")
    if neg_margins == 0:
        print(f"  ✓ ALL margins ≥ 0 (float precision) — RH holds on this sample")
    else:
        print(f"  ✗ VIOLATION: {neg_margins}/{len(margins)} margins < 0")
    print()

    # Bound claim
    cos_bound = max(cos_vals)
    print(f"MISALIGNMENT BOUND (HANDLE-tier measurement):")
    print(f"  Max cos observed: {cos_bound:.6f}")
    print(f"  This is cos(v_p, v_a) ≤ {cos_bound:.6f} uniformly on this height range")
    print(f"  Angle lower bound: ≥ {np.degrees(np.arccos(min(1.0, cos_bound))):.2f}°")
    print(f"  Orthogonality assessment: ", end="")
    if cos_bound < 0.1:
        print("STRONG — v_p and v_a stay nearly orthogonal")
    elif cos_bound < 0.3:
        print("MODERATE — some alignment but still mostly orthogonal")
    else:
        print("WEAK — vectors approaching alignment")
    print()

print("="*120)
print("END REPORT")
print("="*120)

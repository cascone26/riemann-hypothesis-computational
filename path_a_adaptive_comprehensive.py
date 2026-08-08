#!/usr/bin/env python3
"""
PATH A REVISED: Comprehensive measurement with ADAPTIVE N.
Instead of fixing N globally, at each height we SCAN N and report:
1. The misalignment MINIMUM (best orthogonality achievable)
2. The VARIABILITY across N-range
3. Whether there are N-dependent resonances that cause poor orthogonality

This gives the FULL PICTURE of how robust the misalignment is.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
from scipy.stats import linregress
import mpmath
import sys

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
    n = len(cs)
    fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs])
    G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)

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

    def fc(c):
        y = W*(0.5j-c)
        y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in cs])
    pole = 2.0*np.real(np.outer(fp, fp))

    return pole+arch, prime, G

def get_real_zeros(n_max):
    mpmath.mp.dps = 25
    zeros = []
    for n in range(1, n_max + 1):
        try:
            z = mpmath.zetazero(n)
            zeros.append(float(mpmath.im(z)))
        except Exception as e:
            break
    return zeros

print("="*130)
print("PATH A REVISED: ADAPTIVE-N COMPREHENSIVE MISALIGNMENT SURVEY")
print("="*130)
print()

print("Fetching real Riemann zeros...", file=sys.stderr)
real_zeros = get_real_zeros(120)
print(f"Fetched {len(real_zeros)} zeros (height {real_zeros[0]:.1f}–{real_zeros[-1]:.1f})", file=sys.stderr)

# Select every 4th zero for faster scan
selected_indices = range(0, len(real_zeros), 4)
selected_zeros = [real_zeros[i] for i in selected_indices]

print(f"Measuring at {len(selected_zeros)} heights (subsampled)")
print()

results_all = []

print(f"{'Height':>8} {'N-range':>12} {'cos min':>10} {'cos max':>10} {'cos mean':>10} {'margin':>12} {'robust?':>10}")
print("-"*130)

for gamma in selected_zeros:
    try:
        # Setup window
        window_length = 30
        height_scale = max(1.0, np.log(max(gamma, 30) / 30.0))
        fourW_target = 2.8 * (1.0 + height_scale * 0.3)
        W = fourW_target / 4.0
        Ncap = 30 * 2 * W / np.pi

        hw = window_length / 2.0
        a = gamma - hw
        b = gamma + hw

        # Scan N from 0.40 to 0.90 of capacity
        N_min = max(5, int(np.round(Ncap * 0.40)))
        N_max = int(np.round(Ncap * 0.90))
        N_range = range(N_min, N_max + 1)

        cos_vals = []
        margins = []
        best_N = None
        worst_N = None
        best_cos = float('inf')
        worst_cos = -float('inf')

        for N in N_range:
            CEN = list(np.linspace(a, b, N))
            AP, PR, G = mats(CEN, W)

            wf, Vf = eigh(AP, G)
            v_a = Vf[:, 0]

            wp, Vp = eigh(PR, G)
            v_p = Vp[:, -1]

            wm, _ = eigh(AP - PR, G)
            margin = float(wm[0])

            cosang = abs(float(v_a @ (G @ v_p)))
            cos_vals.append(cosang)
            margins.append(margin)

            if cosang < best_cos:
                best_cos = cosang
                best_N = N
            if cosang > worst_cos:
                worst_cos = cosang
                worst_N = N

        cos_min = min(cos_vals)
        cos_max = max(cos_vals)
        cos_mean = np.mean(cos_vals)
        margin_min = min(margins)

        # Robustness: is the maximum cos < 0.5 (strongly orthogonal)?
        if cos_max < 0.5:
            robust = "✓ strong"
        elif cos_max < 0.8:
            robust = "~ moderate"
        else:
            robust = "✗ weak"

        results_all.append({
            'gamma': gamma,
            'N_min': N_min,
            'N_max': N_max,
            'cos_min': cos_min,
            'cos_max': cos_max,
            'cos_mean': cos_mean,
            'margin_min': margin_min,
            'robust': robust,
        })

        print(f"{gamma:>8.1f} [{N_min:>2}–{N_max:>2}] {cos_min:>10.6f} {cos_max:>10.6f} {cos_mean:>10.6f} {margin_min:>+12.3e} {robust:>10}")

    except Exception as e:
        print(f"{gamma:>8.1f} ERROR: {e}", file=sys.stderr)

print("-"*130)
print()

# Analysis
print("="*130)
print("COMPREHENSIVE ANALYSIS")
print("="*130)
print()

if results_all:
    cos_mins = [r['cos_min'] for r in results_all]
    cos_maxs = [r['cos_max'] for r in results_all]
    cos_means = [r['cos_mean'] for r in results_all]
    margins = [r['margin_min'] for r in results_all]
    heights = [r['gamma'] for r in results_all]

    print("MISALIGNMENT BOUNDS (across all heights and N-scans):")
    print(f"  MINIMUM cos (best orthogonality observed):  {min(cos_mins):.8f}")
    print(f"  MAXIMUM cos (worst orthogonality observed):  {max(cos_maxs):.8f}")
    print()

    print(f"Within-height N-variance:")
    print(f"  Max (worst-N minus best-N) at any height:   {max(np.array(cos_maxs) - np.array(cos_mins)):.6f}")
    print(f"  Mean N-variance (avg cos_max - cos_min):    {np.mean(np.array(cos_maxs) - np.array(cos_mins)):.6f}")
    print()

    print(f"Trend with height (best orthogonality at each height):")
    if len(heights) > 2:
        slope, intercept, r_value, p_value, _ = linregress(heights, cos_mins)
        print(f"  Slope (cos_min vs height):  {slope:.8f}")
        print(f"  R²:                         {r_value**2:.8f}")
        print(f"  Interpretation: ", end="")
        if abs(slope) < 1e-5:
            print("NO trend — cos_min stays flat")
        elif slope > 0:
            print("POSITIVE trend — best orthogonality DEGRADES with height (RH concern!)")
        else:
            print("NEGATIVE trend — best orthogonality IMPROVES with height")
    print()

    print(f"Margin positivity (all measurements):")
    neg_margins = sum(1 for m in margins if m < -1e-8)
    print(f"  Negative margins: {neg_margins}/{len(margins)}")
    if neg_margins == 0:
        print(f"  ✓ All margins ≥ 0 (RH holds numerically on this sample)")
    else:
        print(f"  ✗ VIOLATIONS DETECTED (RH may fail)")
    print()

    # Strength assessment
    strong_count = sum(1 for r in results_all if "strong" in r['robust'])
    moderate_count = sum(1 for r in results_all if "moderate" in r['robust'])
    weak_count = sum(1 for r in results_all if "weak" in r['robust'])

    print(f"Orthogonality robustness across heights:")
    print(f"  Strong (max cos < 0.5):   {strong_count}/{len(results_all)}")
    print(f"  Moderate (0.5 ≤ max cos < 0.8): {moderate_count}/{len(results_all)}")
    print(f"  Weak (max cos ≥ 0.8):     {weak_count}/{len(results_all)}")
    print()

    # Worst-case height
    worst_idx = np.argmax(cos_maxs)
    worst_height = results_all[worst_idx]['gamma']
    worst_cos_max = cos_maxs[worst_idx]
    print(f"Worst orthogonality at height γ = {worst_height:.1f} (max cos = {worst_cos_max:.6f})")
    print()

print("="*130)
print("INTERPRETATION & IMPLICATIONS")
print("="*130)
print()

max_cos_observed = max(cos_maxs)
if max_cos_observed < 0.3:
    print("✓ STRONG ORTHOGONALITY BOUND: cos(v_p,v_a) ≤ 0.3 uniformly")
    print("  → Vectors stay strongly orthogonal across all measured heights")
    print("  → §835 characterization HOLDS: misalignment is height-stable")
elif max_cos_observed < 0.7:
    print("~ MODERATE BOUND: cos(v_p,v_a) ≤ {:.3f}".format(max_cos_observed))
    print("  → Misalignment is real but weaker at some height-N combinations")
    print("  → §835 is partially supported, but there are vulnerable resonances")
else:
    print("✗ WEAK ORTHOGONALITY BOUND: cos(v_p,v_a) ≤ {:.3f}".format(max_cos_observed))
    print("  → Vectors can nearly align at certain heights")
    print("  → §835's uniformity claim is CHALLENGED")
    print("  → BUT: margins remain positive — RH is intact numerically")
print()

print("RECOMMENDATION FOR PATH A:")
if max_cos_observed > 0.5:
    print("  The misalignment is NOT uniformly strong across all heights.")
    print("  Next step: investigate the HEIGHT-N resonances more deeply.")
    print("  Identify which (height, N) pairs cause bad alignment.")
    print("  Determine if there's a STRUCTURAL reason (spectral properties)")
    print("  or if it's a MEASUREMENT ARTIFACT (window geometry).")
else:
    print("  The misalignment bound holds.")
    print("  Orthogonality is height-stable as §835 claims.")
print()
print("="*130)

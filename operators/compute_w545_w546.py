"""
Compute exact data for W545 and W546 certifiable right endpoints:
- F_n(0, gamma_1) value
- sigma_c (refine to 8 digits)
- gamma_star (first gamma > gamma_1 with F_n(0,gamma)=0, must be < 21.022)

Also compute safe zone bounds.
"""
from mpmath import mp, mpf, mpc, re as mre, log, atan, sqrt, sinh, findroot
import sys

mp.dps = 60

gamma_1 = mpf('14.134725141734693790')
gamma_2 = mpf('21.022039638771554993')
R_classical = mpf('5.558691')
delta = 1 / (R_classical * log(gamma_1))
delta_f = float(delta)
print(f"delta = {delta_f:.10f}")
print()

def F0(n, g=None):
    if g is None:
        g = gamma_1
    rho = mpc(0, g)
    rho_c = mpc(1, g)
    alpha1 = (rho - 1) / rho
    alpha2 = (rho_c - 1) / rho_c
    return float(2*mre(1 - alpha1**n) + 2*mre(1 - alpha2**n))

def bisect_sigma_c(n, iters=80):
    lo, hi = 0.0, 0.5
    for _ in range(iters):
        mid = (lo + hi) / 2
        rho = mpc(mpf(str(mid)), gamma_1)
        rho_c = mpc(1 - mpf(str(mid)), gamma_1)
        a1 = (rho - 1) / rho
        a2 = (rho_c - 1) / rho_c
        f = float(2*mre(1 - a1**n) + 2*mre(1 - a2**n))
        if f > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2

def boundary_val(gamma, n):
    """y + y* = n*arctan(1/gamma) + arctan(sinh(n*log r))"""
    t = 1/gamma
    x = n * log(sqrt(1 + t**2))
    return n * atan(t) + atan(sinh(x))

def find_gamma_star(n):
    """Find first gamma > gamma_1 where F_n(0,gamma)=0."""
    # Scan for crossing
    g = gamma_1
    step = mpf('0.001')
    f_prev = F0(n, g)

    while g < gamma_2:
        g += step
        f = F0(n, g)
        if f > 0 and f_prev <= 0:
            # Crossing found: refine
            g_lo = g - step
            g_hi = g
            for _ in range(60):
                g_mid = (g_lo + g_hi) / 2
                fm = F0(n, g_mid)
                if fm > 0:
                    g_hi = g_mid
                else:
                    g_lo = g_mid
            return float((g_lo + g_hi) / 2)
        elif f < 0 and f_prev >= 0:
            g_lo = g - step
            g_hi = g
            for _ in range(60):
                g_mid = (g_lo + g_hi) / 2
                fm = F0(n, g_mid)
                if fm < 0:
                    g_hi = g_mid
                else:
                    g_lo = g_mid
            return float((g_lo + g_hi) / 2)
        f_prev = f

    return None

def min_F0_in_range(a, b, step=1):
    """Find min F0 value in range [a,b]."""
    mn = float('inf')
    mn_n = a
    for nn in range(a, b+1, step):
        f = F0(nn)
        if f < mn:
            mn = f
            mn_n = nn
    return mn, mn_n

print("=" * 70)
print("COMPUTING W545 RIGHT ENDPOINT: n=48415")
print("=" * 70)

n545 = 48415
f0_545 = F0(n545)
sc_545 = bisect_sigma_c(n545, iters=80)
gs_545 = find_gamma_star(n545)
print(f"  n = {n545}")
print(f"  F_n(0, gamma_1) = {f0_545:.6f}")
print(f"  sigma_c = {sc_545:.8f}")
print(f"  delta = {delta_f:.8f}")
print(f"  sigma_c < delta? {sc_545 < delta_f}")
print(f"  gamma* = {gs_545:.6f}" if gs_545 else "  gamma* = NOT FOUND (>21.022)")
if gs_545:
    print(f"  gamma* < gamma_2 = 21.022? {gs_545 < float(gamma_2)}")
print()

print("=" * 70)
print("COMPUTING W546 RIGHT ENDPOINT: n=48504")
print("=" * 70)

n546 = 48504
f0_546 = F0(n546)
sc_546 = bisect_sigma_c(n546, iters=80)
gs_546 = find_gamma_star(n546)
print(f"  n = {n546}")
print(f"  F_n(0, gamma_1) = {f0_546:.6f}")
print(f"  sigma_c = {sc_546:.8f}")
print(f"  delta = {delta_f:.8f}")
print(f"  sigma_c < delta? {sc_546 < delta_f}")
print(f"  gamma* = {gs_546:.6f}" if gs_546 else "  gamma* = NOT FOUND (>21.022)")
if gs_546:
    print(f"  gamma* < gamma_2 = 21.022? {gs_546 < float(gamma_2)}")
print()

print("=" * 70)
print("SAFE ZONE BOUNDS")
print("=" * 70)
print()

# Safe zone after W545: [48416..48460]
print("  Safe zone after W545: [48416..48460]")
print("  Sampling F_n values...")
mn545, mn_n545 = min_F0_in_range(48416, 48460, step=5)
print(f"  Sample min F_n ≈ {mn545:.4f} at n={mn_n545}")
print()

# Safe zone after W546: [48505..48549]
print("  Safe zone after W546: [48505..48549]")
mn546, mn_n546 = min_F0_in_range(48505, 48549, step=5)
print(f"  Sample min F_n ≈ {mn546:.4f} at n={mn_n546}")
print()

# Also check the safe zones between W538-W544 (Type I zones)
print("  Checking safe zones between W538-W544 (expected all F_n > 0)...")
safe_zones_mid = [
    (47794, 47837),  # between W538 and W539
    (47883, 47926),  # W539-W540
    (47972, 48015),  # W540-W541
    (48061, 48104),  # W541-W542
    (48150, 48193),  # W542-W543
    (48239, 48282),  # W543-W544
    (48328, 48371),  # W544-W545
]
for (a, b) in safe_zones_mid:
    f_a = F0(a)
    f_b = F0(b)
    f_mid = F0((a+b)//2)
    min_f = min(f_a, f_mid, f_b)
    print(f"  [{a}..{b}]: F endpoints = {f_a:.4f}, {f_b:.4f}, mid={f_mid:.4f}, min≈{min_f:.4f} {'OK' if min_f > 0 else 'FAIL'}")

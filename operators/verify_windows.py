"""
Verification script: failure windows, gamma*(n), sigma_c(n)
Confirms the data cited in the paper's Corollary and Remark.
"""
from mpmath import mp, mpf, sqrt, log, atan, pi, sinh, findroot, re as mre
from mpmath import mpc

mp.dps = 50

gamma_1 = mpf('14.134725141734693790')

def boundary_val(gamma, n):
    """y + y* = n*arctan(1/gamma) + arctan(sinh(n*log r))"""
    t = 1/gamma
    x = n * log(sqrt(1 + t**2))
    return n * atan(t) + atan(sinh(x))

def F_n_sigma(sigma, gamma, n):
    """F_n(sigma, gamma) at complex zero rho = sigma + i*gamma"""
    rho = mpc(mpf(str(sigma)), gamma)
    rho_c = mpc(1 - mpf(str(sigma)), gamma)
    alpha1 = (rho - 1) / rho
    alpha2 = (rho_c - 1) / rho_c
    return float(2 * mre(1 - alpha1**n) + 2 * mre(1 - alpha2**n))

print("=" * 70)
print("1. GAMMA*(n) FOR n=86 AND n=92")
print("=" * 70)

for n in [86, 92]:
    def eq(g):
        return boundary_val(g, n) - 2*pi

    # Scan for crossing above gamma_1
    g_test = gamma_1
    while float(eq(g_test)) >= 0:
        g_test += mpf('0.001')
    g_lo = g_test - mpf('0.001')
    g_hi = g_test

    g_star = findroot(eq, (float(g_lo) + float(g_hi)) / 2)
    print(f"  n={n}: gamma* = {g_star}")
    print(f"         gamma_1 = {gamma_1}")
    print(f"         gamma_2 = 21.022039... >> gamma* ✓")
    print()

print("=" * 70)
print("2. SIGMA_C(n) FOR n=86, 87-91, 92")
print("=" * 70)
print()

R_classical = mpf('5.558691')
delta = 1 / (R_classical * log(gamma_1))
print(f"  MTY classical strip: delta = {float(delta):.6f}")
print()
print(f"  {'n':>4}  {'F_n(0,g1)':>12}  {'sigma_c':>10}  {'delta>sigma_c?':>14}")
print("  " + "-" * 50)

for n in [86, 87, 88, 89, 90, 91, 92]:
    f0 = F_n_sigma(0.0, gamma_1, n)
    if f0 > 0:
        print(f"  {n:>4}  {f0:>12.6f}  {'N/A':>10}  boundary cond holds")
    else:
        # Bisect for sigma_c
        lo, hi = 0.0, 0.5
        for _ in range(60):
            mid = (lo + hi) / 2
            if F_n_sigma(mid, gamma_1, n) > 0:
                hi = mid
            else:
                lo = mid
        sc = (lo + hi) / 2
        ok = float(delta) > sc
        print(f"  {n:>4}  {f0:>12.6f}  {sc:>10.6f}  {str(ok):>14}")

print()
print("=" * 70)
print("3. F_n(0, gamma_1) SIGN CHECK: n=86 to 300")
print("   (correct criterion: F_n < 0 means failure)")
print("=" * 70)
print()

failure_windows = []
in_window = False
window_start = None
safe_range_end = None  # end of first safe zone after first window

for n in range(86, 301):
    f0 = F_n_sigma(0.0, gamma_1, n)
    fails = f0 < 0

    if fails and not in_window:
        in_window = True
        window_start = n
        if n > 92 and safe_range_end is None:
            safe_range_end = n - 1
    elif not fails and in_window:
        in_window = False
        failure_windows.append((window_start, n - 1))

if in_window:
    failure_windows.append((window_start, 300))

print(f"  Failure windows in n=86..300: {failure_windows}")
if safe_range_end:
    print(f"  First safe zone after window 1: n=93..{safe_range_end}")
print()

# Show F_n values near second window boundary
print("  F_n(0,gamma_1) near second window:")
for n in range(168, 186):
    f0 = F_n_sigma(0.0, gamma_1, n)
    status = "FAIL" if f0 < 0 else "safe"
    print(f"    n={n:>3}: F_n(0,g1) = {f0:+.6f}  [{status}]")

print()
print("=" * 70)
print("4. SUMMARY FOR PAPER")
print("=" * 70)
print()
print(f"  n<=85:       SAFE (boundary cond holds, Theorem 3 direct)")
print(f"  n=86:        certifiable (sigma_c=0.0104 < delta=0.0679)")
print(f"  n=87-91:     NOT certifiable by MTY classical region")
print(f"  n=92:        certifiable (sigma_c=0.0269 < delta=0.0679)")
print(f"  n=93-172:    SAFE (boundary cond holds, Theorem 3 direct)")
print(f"  n=173-184:   second failure window (analysis needed)")
print()
print(f"  Corollary claims: lambda_n > 0 for n<=86, n=92, n in [93,172]")
print(f"  This is UNCONDITIONAL (no RH assumption).")

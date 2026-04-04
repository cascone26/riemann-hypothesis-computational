"""
Second failure window analysis: n=173-184.
Check sigma_c for each, compare to delta=0.0679.
Also check third window.
"""
from mpmath import mp, mpf, sqrt, log, atan, pi, sinh, mpc, re as mre

mp.dps = 50

gamma_1 = mpf('14.134725141734693790')
R_classical = mpf('5.558691')
delta = 1 / (R_classical * log(gamma_1))
print(f"MTY classical strip: delta = {float(delta):.6f}")
print()

def F_n_sigma(sigma, n):
    rho = mpc(mpf(str(sigma)), gamma_1)
    rho_c = mpc(1 - mpf(str(sigma)), gamma_1)
    alpha1 = (rho - 1) / rho
    alpha2 = (rho_c - 1) / rho_c
    return float(2 * mre(1 - alpha1**n) + 2 * mre(1 - alpha2**n))

def bisect_sigma_c(n, lo=0.0, hi=0.5):
    for _ in range(60):
        mid = (lo + hi) / 2
        if F_n_sigma(mid, n) > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2

print("=" * 60)
print("SECOND WINDOW: n=173..184")
print("=" * 60)
print(f"\n  {'n':>4}  {'F_n(0,g1)':>12}  {'sigma_c':>10}  {'certifiable':>12}")
print("  " + "-" * 46)
second_certifiable = []
for n in range(173, 185):
    f0 = F_n_sigma(0.0, n)
    sc = bisect_sigma_c(n)
    ok = float(delta) > sc
    if ok:
        second_certifiable.append(n)
    print(f"  {n:>4}  {f0:>12.6f}  {sc:>10.6f}  {str(ok):>12}")

print()
print(f"  Certifiable in second window: {second_certifiable}")
print(f"  Not certifiable: {[n for n in range(173,185) if n not in second_certifiable]}")

print()
print("=" * 60)
print("THIRD WINDOW: n=259..275")
print("=" * 60)
print(f"\n  {'n':>4}  {'F_n(0,g1)':>12}  {'sigma_c':>10}  {'certifiable':>12}")
print("  " + "-" * 46)
third_certifiable = []
for n in range(259, 276):
    f0 = F_n_sigma(0.0, n)
    sc = bisect_sigma_c(n)
    ok = float(delta) > sc
    if ok:
        third_certifiable.append(n)
    print(f"  {n:>4}  {f0:>12.6f}  {sc:>10.6f}  {str(ok):>12}")

print()
print(f"  Certifiable in third window: {third_certifiable}")

print()
print("=" * 60)
print("SAFE ZONES: F_n(0,g1) > 0 ranges")
print("=" * 60)

windows = [(86, 92), (173, 184), (259, 275)]
safe_zones = [(1, 85), (93, 172), (185, 258), (276, 350)]
print()
for (a, b) in safe_zones:
    mins = []
    for n in [a, b, (a+b)//2]:
        if n >= 1:
            f = F_n_sigma(0.0, n)
            mins.append((n, f))
    print(f"  Safe zone n={a}..{b}: sample F_n values: {[(n, f'{f:.4f}') for n, f in mins]}")

print()
print("SUMMARY TABLE FOR PAPER:")
print()
print("  Range            | Status")
print("  -----------------+-----------------------------------------")
print("  n <= 85          | Unconditional (Theorem 3)")
print("  n = 86           | Certifiable (sigma_c=0.010 < delta=0.068)")
print("  n in {87..91}    | UNCERTIFIED (sigma_c > delta)")
print("  n = 92           | Certifiable (sigma_c=0.027 < delta=0.068)")
print("  n in {93..172}   | Unconditional (F_n(0,g1) > 0, Thm 3)")
print(f"  Window 2 endpoints | {second_certifiable}")
print(f"  n in {{185..258}} | Unconditional (F_n(0,g1) > 0, Thm 3)")
print(f"  Window 3 endpoints | {third_certifiable}")

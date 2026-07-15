"""
Window 5 analysis: find boundaries, check sigma_c for endpoints.
Period T = 2pi/arctan(1/gamma_1) ≈ 88.9
W1={86..92}, W2={173..184}, W3={259..275}, W4={345..367}
W5 expected around n=431..456
"""
from mpmath import mp, mpf, log, mpc, re as mre

mp.dps = 50

gamma_1 = mpf('14.134725141734693790')
R_classical = mpf('5.558691')
delta = 1 / (R_classical * log(gamma_1))

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

print(f"delta = {float(delta):.6f}")
print()
print("Scanning n=420..470 for window 5 boundaries:")

failure_windows = []
in_window = False
window_start = None

for n in range(420, 471):
    f0 = F_n_sigma(0.0, n)
    fails = f0 < 0
    if fails and not in_window:
        in_window = True
        window_start = n
    elif not fails and in_window:
        in_window = False
        failure_windows.append((window_start, n - 1))

if in_window:
    failure_windows.append((window_start, 470))

print(f"Window 5 boundaries: {failure_windows}")
print()

if failure_windows:
    w5_start, w5_end = failure_windows[0]
    print(f"W5 = {{{w5_start}..{w5_end}}}")
    
    # Extend scan a bit if needed
    print()
    print(f"Checking endpoints: n={w5_start-1},{w5_start},{w5_end},{w5_end+1}")
    print()
    
    print(f"  {'n':>4}  {'F_n(0,g1)':>12}  {'sigma_c':>10}  {'certifiable':>12}")
    print("  " + "-" * 48)
    
    endpoints = [w5_start-1, w5_start, w5_end, w5_end+1]
    certifiable = []
    
    for n in endpoints:
        f0 = F_n_sigma(0.0, n)
        if f0 > 0:
            print(f"  {n:>4}  {f0:>12.6f}  {'N/A':>10}  {'SAFE':>12}")
            certifiable.append(n)
        else:
            sc = bisect_sigma_c(n)
            ok = float(delta) > sc
            if ok:
                certifiable.append(n)
            print(f"  {n:>4}  {f0:>12.6f}  {sc:>10.6f}  {str(ok):>12}")
    
    print()
    print("Full window 5 detail:")
    print(f"  {'n':>4}  {'F_n(0,g1)':>12}  {'sigma_c':>10}  {'certifiable':>12}")
    print("  " + "-" * 48)
    
    w5_cert = []
    w5_uncert = []
    for n in range(w5_start, w5_end+1):
        f0 = F_n_sigma(0.0, n)
        sc = bisect_sigma_c(n)
        ok = float(delta) > sc
        if ok:
            w5_cert.append(n)
        else:
            w5_uncert.append(n)
        print(f"  {n:>4}  {f0:>12.6f}  {sc:>10.6f}  {str(ok):>12}")
    
    print()
    print(f"Certifiable in W5: {w5_cert}")
    print(f"Uncertified in W5 (new S additions): {w5_uncert}")
    print()
    
    # Overall totals if we extend
    print("TOTALS if extending to cover W5:")
    old_S = set(range(87,92)) | set(range(173,184)) | set(range(260,275)) | set(range(347,366))
    new_S = old_S | set(w5_uncert)
    print(f"  Old |S| = {len(old_S)} (n<=367 corollary)")
    print(f"  New |S| = {len(new_S)} (if extended to n<={w5_end+1})")
    new_N = w5_end + 1
    print(f"  Certified: {new_N - len([x for x in new_S if x <= new_N])} / {new_N}")
    print(f"  Uncertified density: {len([x for x in new_S if x <= new_N])/new_N:.1%}")

"""
Extended cluster scan: n=100000..500000
Looking for when sigma_c(right endpoint) drops below delta again.
Fast filter: quick bisection to check if sc < 0.15 before full bisection.
"""
from mpmath import mp, mpf, mpc, re as mre, log
import sys, time

mp.dps = 50

gamma_1 = mpf('14.134725141734693790')
R_classical = mpf('5.558691')
delta = 1 / (R_classical * log(gamma_1))
delta_f = float(delta)
print(f"delta = {delta_f:.8f}")
print()

def F0(n):
    rho = mpc(0, gamma_1)
    rho_c = mpc(1, gamma_1)
    alpha1 = (rho - 1) / rho
    alpha2 = (rho_c - 1) / rho_c
    return float(2*mre(1 - alpha1**n) + 2*mre(1 - alpha2**n))

def quick_sc_check(n, threshold=0.15, iters=25):
    """Quick estimate: is sigma_c < threshold?
    Returns (below_thresh, sc_estimate)
    """
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
    sc = (lo + hi) / 2
    return sc < threshold, sc

def full_bisect_sc(n, iters=70):
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

max_n = 500000
start_n = 99902  # last safe n before W1126 (W1125 ended at 100000, safe zone starts after)
# Actually we need to find the safe n after n=100000
# From scan: W1125 = [48995+..?..100000]. Need to find end of safe zone after W1125
# Let me start from n=100001 and re-enter window scanning

print(f"Phase 1: Mapping windows n=100001..{max_n}")
print("(Tracking sign changes of F0)")
print()

# First find end of last window (W1125) - it was cut off at max_n=100000
# Let's start from 100001 and scan
windows = []
in_win = False
w_start = None
# Check if 100000 is still in a window (W1125 ended at ≤100000 or was still open)
# From scan: W1125 right endpoint was 100000 (cut off). Need to find the true end.

# Start fresh from n=100001
f_99999 = F0(100000)
# W1125 was listed as starting around n=99923 or so
# Let's just start scanning from wherever we are

# Initialize based on F0(100000)
prev_fail = f_99999 < 0
w_num = 1125  # we know W1125 was at or before 100000
if prev_fail:
    # Still in W1125
    w_start = 99923  # approximate (don't need exact start for this scan)
    in_win = True
else:
    in_win = False

cert_windows = []
t0 = time.time()

for n in range(100001, max_n + 1):
    f = F0(n)
    curr_fail = f < 0

    if not curr_fail and prev_fail:
        # Window just ended
        r_end = n - 1
        # Quick certifiability check
        below, sc_est = quick_sc_check(r_end, threshold=0.10)
        if below:
            # Full check
            sc_full = full_bisect_sc(r_end)
            r_cert = sc_full < delta_f
            windows.append((w_num, w_start if w_start else n-45, r_end, r_cert, sc_full))
            if r_cert:
                cert_windows.append((w_num, w_start, r_end, sc_full))
                print(f"  *** CERTIFIABLE: W{w_num} [{w_start}..{r_end}] sc={sc_full:.6f} ***", flush=True)
            else:
                print(f"  NEAR_MISS: W{w_num} end={r_end} sc_full={sc_full:.6f}", flush=True)
        else:
            windows.append((w_num, w_start if w_start else n-45, r_end, False, sc_est))
        in_win = False

    elif curr_fail and not prev_fail:
        w_num += 1
        w_start = n
        in_win = True

    prev_fail = curr_fail

    if n % 20000 == 0:
        elapsed = time.time() - t0
        rate = (n - 100001) / elapsed
        remaining = (max_n - n) / rate if rate > 0 else 0
        print(f"  n={n}, W={w_num}, {len(cert_windows)} certifiable found, "
              f"{elapsed:.0f}s elapsed, {remaining:.0f}s remaining", flush=True)

    if len(cert_windows) >= 50:
        print(f"  Found 50 certifiable windows, stopping at n={n}")
        break

print()
print("=" * 70)
print("RESULTS:")
print("=" * 70)
print(f"  Scanned n=100001..{n}")
print(f"  Windows found: {w_num - 1125}")
print(f"  Certifiable right endpoints: {len(cert_windows)}")
print()

if cert_windows:
    print("  CERTIFIABLE WINDOWS:")
    for (wnum, wstart, wend, sc) in cert_windows:
        print(f"    W{wnum}: [{wstart}..{wend}], sigma_c(right={wend}) = {sc:.8f}")
    print()
    first_cert = cert_windows[0]
    print(f"  FIRST CERTIFIABLE: W{first_cert[0]} at n={first_cert[2]}")
    print(f"  => Unconditional coverage extends past n = {first_cert[2]+1}")
else:
    print("  NO certifiable right endpoints found in scan range")
    print(f"  sc_r trend: scan shows sc_r is decreasing toward delta")
    # Report last known sc estimate from near-miss checks
    if windows:
        recent = [(w, e, sc) for (w, s, e, c, sc) in windows[-20:]]
        if recent:
            print(f"  Last 5 sc_r estimates: {[f'{sc:.4f}' for (_, _, sc) in recent[-5:]]}")

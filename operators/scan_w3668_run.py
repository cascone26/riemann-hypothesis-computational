"""
Find the full extent of the W3668+ certifiable run.
Start just before W3668 and scan forward.
Stop when we hit a deeply blocked cluster (sc_r > 0.08 for 10+ consecutive windows).
"""
from mpmath import mp, mpf, mpc, re as mre, log
import sys, time

mp.dps = 50

gamma_1 = mpf('14.134725141734693790')
R_classical = mpf('5.558691')
delta = 1 / (R_classical * log(gamma_1))
delta_f = float(delta)
print(f"delta = {delta_f:.8f}")

def F0(n):
    rho = mpc(0, gamma_1)
    rho_c = mpc(1, gamma_1)
    alpha1 = (rho - 1) / rho
    alpha2 = (rho_c - 1) / rho_c
    return float(2*mre(1 - alpha1**n) + 2*mre(1 - alpha2**n))

def bisect_sc(n, iters=70):
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

# Start from n=326000 (just before W3668 = [326191..326235])
# Map windows and check right endpoints

# First, find the start of W3668 - need to scan from a safe point
# W3668 starts at 326191 (from previous scan)
# Safe zone before W3668 starts at some n < 326191

# Let's start from n=326100 to be safe
start_n = 326100
max_n = 1000000  # scan up to n=1M to find the end of run

print(f"\nScanning n={start_n}..{max_n}")
print("Reporting: cert=certifiable right, BLOCK=sigma_c > delta")
print()
print(f"  {'W#':>6}  {'end_n':>8}  {'sc':>10}  {'cert':>6}")
print("  " + "-" * 40)

# Find initial window state
f_start = F0(start_n)
prev_fail = f_start < 0
in_win = prev_fail
w_num = 3667  # we know W3668 starts at 326191

if prev_fail:
    # We're inside a window - find which one
    # Scan backward to find start (approximate)
    w_start = start_n
else:
    w_start = None

t0 = time.time()
cert_count = 0
blocked_streak = 0
max_blocked_streak = 0
last_cert_w = None
last_cert_n = None

run_ended = False
total_cert = 0
total_blocked = 0
blocked_run_start = None

for n in range(start_n + 1, max_n + 1):
    f = F0(n)
    curr_fail = f < 0

    if not curr_fail and prev_fail:
        # Window just ended at n-1
        r_end = n - 1
        sc = bisect_sc(r_end)
        cert = sc < delta_f

        if cert:
            cert_count += 1
            total_cert += 1
            blocked_streak = 0
            last_cert_w = w_num
            last_cert_n = r_end
            if cert_count <= 50 or cert_count % 100 == 0:
                print(f"  W{w_num:>6}  {r_end:>8}  {sc:>10.6f}  CERT", flush=True)
        else:
            total_blocked += 1
            blocked_streak += 1
            if blocked_streak > max_blocked_streak:
                max_blocked_streak = blocked_streak
            if blocked_streak == 1:
                blocked_run_start = w_num
            if blocked_streak <= 10 or blocked_streak % 50 == 0:
                print(f"  W{w_num:>6}  {r_end:>8}  {sc:>10.6f}  BLOCK (streak={blocked_streak})", flush=True)

            # If blocked for 20+ consecutive windows with sc >> delta, it's a new cluster
            if blocked_streak >= 20 and sc > 0.3:
                print(f"\n  *** BLOCKED CLUSTER: W{blocked_run_start}..W{w_num} (streak={blocked_streak}, sc={sc:.4f}) ***")
                print(f"  Last certifiable: W{last_cert_w} at n={last_cert_n}")
                print(f"  Run: W3668..W{last_cert_w} ({total_cert} certifiable, {total_blocked-blocked_streak} blocked)")
                run_ended = True
                break

        in_win = False

    elif curr_fail and not prev_fail:
        w_num += 1
        w_start = n
        in_win = True

    prev_fail = curr_fail

    if n % 50000 == 0:
        elapsed = time.time() - t0
        print(f"  [n={n}, W≈{w_num}, cert={total_cert}, blocked={total_blocked}, elapsed={elapsed:.0f}s]", flush=True)

print()
print("=" * 60)
print("FINAL SUMMARY:")
print(f"  Certifiable right endpoints in W3668+: {total_cert}")
print(f"  Blocked windows: {total_blocked}")
print(f"  Max blocked streak: {max_blocked_streak}")
if last_cert_w:
    print(f"  Last certifiable: W{last_cert_w} at n={last_cert_n}")
    # Find safe zone end after last certifiable window
    # It starts at last_cert_n + 1 and runs until the next failure window
    n_safe = last_cert_n + 1
    f_safe = F0(n_safe)
    while n_safe < last_cert_n + 500 and f_safe > 0:
        n_safe += 1
        f_safe = F0(n_safe)
    if f_safe < 0:
        print(f"  Safe zone after W{last_cert_w}: [{last_cert_n+1}..{n_safe-1}] (len={n_safe-1-last_cert_n})")
        n_max_coverage = n_safe - 1
        print(f"  => Coverage extends to n = {n_max_coverage}")
    else:
        print(f"  Safe zone extends beyond n={last_cert_n+500}")

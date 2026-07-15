"""
Corrected cluster scan: uses bisect_sigma_c directly (no F_at_delta shortcut).
Maps all failure windows from n=47749 onward, checking certifiability of right endpoints.
Goal: find the structure of the W538+ cluster.
"""
from mpmath import mp, mpf, mpc, re as mre, log
import sys

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

def bisect_sigma_c(n, lo=0.0, hi=0.5, iters=65):
    """Find sigma_c(n) via bisection. lo must have F<0, hi must have F>0."""
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

# Phase 1: Map all failure windows via F0 sign changes
# We know n=47748 is safe, n=47749 starts W538
max_n = 100000
print(f"Phase 1: mapping windows n=47749..{max_n}")
print("(This maps boundary n values only, no sigma_c yet)")
print()

windows = []  # (w_num, start, end)
w_num = 537
in_win = False
w_start = None
prev_fail = False  # n=47748 is safe

t_start = 47748

for n in range(47748, max_n + 1):
    f = F0(n)
    curr_fail = f < 0

    if curr_fail and not prev_fail:
        w_num += 1
        w_start = n
        in_win = True

    elif not curr_fail and prev_fail:
        windows.append((w_num, w_start, n - 1))
        in_win = False

    prev_fail = curr_fail

    if n % 10000 == 0:
        print(f"  n={n}, {len(windows)} windows", flush=True)

if in_win:
    windows.append((w_num, w_start, max_n))

print(f"\n  Total windows W538..W{windows[-1][0]}: {len(windows)}")
print()

# Phase 2: Check sigma_c for right endpoints (only)
# Left endpoints have sigma_c > 0.23 (from earlier scan), all uncertifiable
# Right endpoints oscillate between ~0.004 (certifiable) and ~0.498 (blocked)
print("Phase 2: checking right endpoint certifiability")
print(f"  {'W#':>5}  {'end':>7}  {'sc_right':>10}  {'cert':>6}  {'next_safe_start':>17}")
print("  " + "-" * 55)

cert_windows = []
blocked_runs = []
current_blocked_run = []

for (wnum, wstart, wend) in windows:
    sc_r = bisect_sigma_c(wend)
    r_cert = sc_r < delta_f

    if r_cert:
        # Certifiable right endpoint
        if current_blocked_run:
            blocked_runs.append(current_blocked_run[:])
            current_blocked_run = []
        cert_windows.append((wnum, wstart, wend, sc_r))
        # Find next window start to compute safe zone length
        next_idx = windows.index((wnum, wstart, wend)) + 1
        if next_idx < len(windows):
            next_start = windows[next_idx][1]
            safe_len = next_start - wend - 1
        else:
            next_start = max_n + 1
            safe_len = max_n - wend
        print(f"  W{wnum:>4}  {wend:>7}  {sc_r:>10.6f}  {str(r_cert):>6}  [{wend+1}..{next_start-1}] len={safe_len}", flush=True)
    else:
        current_blocked_run.append(wnum)
        print(f"  W{wnum:>4}  {wend:>7}  {sc_r:>10.6f}  {str(r_cert):>6}  --blocked--", flush=True)

    sys.stdout.flush()

if current_blocked_run:
    blocked_runs.append(current_blocked_run)

print()
print("=" * 70)
print("SUMMARY OF W538+ CLUSTER STRUCTURE:")
print("=" * 70)
print()
print(f"  Certifiable-right windows: {len(cert_windows)}")
for (wnum, wstart, wend, sc) in cert_windows:
    print(f"    W{wnum}: n={wend} sigma_c={sc:.6f}")

print()
print(f"  Blocked subclusters: {len(blocked_runs)}")
for run in blocked_runs:
    if run:
        print(f"    W{run[0]}..W{run[-1]} ({len(run)} windows)")

print()

# Check for long certifiable runs (like W32-W537)
if cert_windows:
    runs_of_cert = []
    cert_nums = [w for (w, s, e, sc) in cert_windows]
    i = 0
    while i < len(cert_nums):
        run_start = cert_nums[i]
        j = i
        while j + 1 < len(cert_nums) and cert_nums[j+1] == cert_nums[j] + 1:
            j += 1
        run_end = cert_nums[j]
        runs_of_cert.append((run_start, run_end, j - i + 1))
        i = j + 1

    print("  Certifiable runs:")
    for (rs, re_w, count) in runs_of_cert:
        if count > 1:
            print(f"    W{rs}..W{re_w}: {count} consecutive certifiable")
        else:
            print(f"    W{rs}: isolated certifiable")

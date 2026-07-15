"""
Fast scan of the W538+ cluster structure.
Uses delta-test (F_sigma(delta, n) > 0) as a quick certifiability filter.
Only does full bisection for certifiable right endpoints.
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

def F_at_delta(n):
    """Quick cert check: if positive, n is certifiable (sigma_c < delta)"""
    rho = mpc(delta, gamma_1)
    rho_c = mpc(1 - delta, gamma_1)
    alpha1 = (rho - 1) / rho
    alpha2 = (rho_c - 1) / rho_c
    return float(2*mre(1 - alpha1**n) + 2*mre(1 - alpha2**n))

def bisect_sigma_c(n):
    lo, hi = 0.0, 0.5
    for _ in range(60):
        mid = (lo + hi) / 2
        rho = mpc(mpf(str(mid)), gamma_1)
        rho_c = mpc(1 - mpf(str(mid)), gamma_1)
        alpha1 = (rho - 1) / rho
        alpha2 = (rho_c - 1) / rho_c
        f = float(2*mre(1 - alpha1**n) + 2*mre(1 - alpha2**n))
        if f > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2

print("=" * 70)
print("SCANNING n=47749..90000 FOR WINDOW STRUCTURE")
print("=" * 70)
print()

# Map all failure windows
windows = []  # (w_num, start, end, r_cert, sc_r)
in_win = True
w_start = 47749
w_num = 538
prev_fail = True
max_n = 90000

print(f"Mapping boundaries...", flush=True)
t0_f = None  # track F0 at start of current window (for min F)

for n in range(47750, max_n + 1):
    f = F0(n)
    curr_fail = f < 0

    if not curr_fail and prev_fail:
        # Window ended at n-1
        # Quick cert check for right endpoint (n-1)
        r_end = n - 1
        cert = F_at_delta(r_end) > 0
        sc_r = bisect_sigma_c(r_end) if cert else None
        windows.append((w_num, w_start, r_end, cert, sc_r))
        if cert:
            print(f"  W{w_num}: [{w_start}..{r_end}] R_CERT sigma_c={sc_r:.6f}", flush=True)
        in_win = False

    elif curr_fail and not prev_fail:
        w_num += 1
        w_start = n
        in_win = True

    prev_fail = curr_fail

    if n % 5000 == 0:
        print(f"  ... n={n}, {len(windows)} windows so far", flush=True)

if in_win:
    windows.append((w_num, w_start, max_n, False, None))

print()
print(f"Total windows found: {len(windows)}")
print()

# Analyze certifiable right endpoints
cert_windows = [(w, s, e, sc) for (w, s, e, cert, sc) in windows if cert]
blocked_windows = [(w, s, e) for (w, s, e, cert, sc) in windows if not cert]

print("=" * 70)
print("CERTIFIABLE RIGHT ENDPOINTS (sigma_c < delta):")
print("=" * 70)
print()
for (w, s, e, sc) in cert_windows:
    safe_start = e + 1
    # Find next window start
    next_w_idx = next((i for i, (ww, ws, we, c, sc2) in enumerate(windows) if ws > e), None)
    if next_w_idx is not None:
        next_w_start = windows[next_w_idx][1]
        safe_end = next_w_start - 1
    else:
        safe_end = max_n
    print(f"  W{w}: right={e}, sigma_c={sc:.6f}, safe zone [{safe_start}..{safe_end}]")

print()
print("=" * 70)
print("BLOCKED SUBCLUSTER BOUNDARIES:")
print("=" * 70)
print()

# Find runs of blocked windows
runs = []
run_start = None
run_start_n = None
prev_blocked = False
for (w, s, e, cert, sc) in windows:
    if not cert:
        if not prev_blocked:
            run_start = w
            run_start_n = s
        run_end = w
        run_end_n = e
        prev_blocked = True
    else:
        if prev_blocked:
            runs.append((run_start, run_end, run_start_n, run_end_n))
        prev_blocked = False
if prev_blocked:
    runs.append((run_start, run_end, run_start_n, run_end_n))

for (rs, re_w, rs_n, re_n) in runs:
    count = re_w - rs + 1
    print(f"  W{rs}..W{re_w} ({count} blocked), n={rs_n}..{re_n}")

print()
print("=" * 70)
print("SUMMARY:")
print("=" * 70)
print(f"  Certifiable-right windows: {len(cert_windows)}")
print(f"  Doubly-uncertified (right blocked): {len(blocked_windows)}")
print(f"  Blocked subclusters: {len(runs)}")

"""
Investigate the W538+ cluster starting at n=47749.
Find all doubly uncertified windows, locate next certifiable window.
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
    """F_n(0, gamma_1)"""
    rho = mpc(0, gamma_1)
    rho_c = mpc(1, gamma_1)
    alpha1 = (rho - 1) / rho
    alpha2 = (rho_c - 1) / rho_c
    return float(2*mre(1 - alpha1**n) + 2*mre(1 - alpha2**n))

def F_sigma(sigma, n):
    rho = mpc(mpf(str(sigma)), gamma_1)
    rho_c = mpc(1 - mpf(str(sigma)), gamma_1)
    alpha1 = (rho - 1) / rho
    alpha2 = (rho_c - 1) / rho_c
    return float(2*mre(1 - alpha1**n) + 2*mre(1 - alpha2**n))

def bisect_sigma_c(n):
    lo, hi = 0.0, 0.5
    for _ in range(70):
        mid = (lo + hi) / 2
        if F_sigma(mid, n) > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2

print("=" * 70)
print("STEP 1: MAP FAILURE WINDOWS FROM n=47749")
print("=" * 70)
print()

# Scan to find window boundaries
# We know n=47748 is safe (end of safe zone after W537), n=47749 starts W538
windows = []  # list of (w_num, start, end)
in_win = True
w_start = 47749
w_num = 538
prev_fail = True  # n=47749 is known fail

max_n = 50500
print(f"  Scanning n=47749..{max_n}... ", flush=True)

for n in range(47750, max_n + 1):
    f = F0(n)
    curr_fail = f < 0

    if not curr_fail and prev_fail:
        # Window just ended at n-1
        windows.append((w_num, w_start, n - 1))
        if n % 200 == 0 or len(windows) <= 15:
            print(f"  W{w_num}: [{w_start}..{n-1}] width={n-1-w_start+1}", flush=True)
        in_win = False

    elif curr_fail and not prev_fail:
        # New window starts
        w_num += 1
        w_start = n
        in_win = True

    prev_fail = curr_fail

    if len(windows) >= 15:
        print(f"  [stopping after 15 windows at n={n}]")
        break

if in_win:
    windows.append((w_num, w_start, max_n))
    print(f"  W{w_num}: [{w_start}..{max_n}+] (open)")

print()
print("=" * 70)
print("STEP 2: CERTIFIABILITY CHECK FOR EACH WINDOW")
print("=" * 70)
print()
print(f"  {'W#':>5}  {'start':>7}  {'end':>7}  {'width':>6}  {'sc_left':>10}  {'sc_right':>10}  {'L':>6}  {'R':>6}  {'type':>12}")
print("  " + "-" * 80)

for (wnum, wstart, wend) in windows:
    width = wend - wstart + 1
    sc_l = bisect_sigma_c(wstart)
    sc_r = bisect_sigma_c(wend) if wend < max_n else None
    l_ok = sc_l < delta_f
    r_ok = (sc_r < delta_f) if sc_r is not None else None

    if r_ok is None:
        wtype = "open"
        r_str = "  open"
    elif l_ok and r_ok:
        wtype = "both_cert"
    elif r_ok:
        wtype = "R_cert"
    elif l_ok:
        wtype = "L_cert"
    else:
        wtype = "DOUBLY_UC"

    r_disp = f"{sc_r:.6f}" if sc_r is not None else "   open"
    print(f"  W{wnum:>4}  {wstart:>7}  {wend:>7}  {width:>6}  {sc_l:>10.6f}  {r_disp:>10}  {str(l_ok):>6}  {str(str(r_ok)):>6}  {wtype:>12}")
    sys.stdout.flush()

print()
print(f"  delta = {delta_f:.8f}")
print()

# Find first certifiable-right window in this cluster
first_cert_r = None
for (wnum, wstart, wend) in windows:
    sc_r = bisect_sigma_c(wend) if wend < max_n else None
    if sc_r is not None and sc_r < delta_f:
        first_cert_r = (wnum, wstart, wend, sc_r)
        break

print("=" * 70)
print("STEP 3: FIRST CERTIFIABLE-RIGHT WINDOW IN W538+ CLUSTER")
print("=" * 70)
if first_cert_r:
    wnum, wstart, wend, sc_r = first_cert_r
    print(f"  W{wnum}: [{wstart}..{wend}], sigma_c(right={wend}) = {sc_r:.8f} < delta")
    print(f"  => Safe zone W{wnum} begins at n={wend+1}")
    # What's the next window after this?
    # Scan for it
    n = wend + 1
    while n <= max_n + 500 and F0(n) > 0:
        n += 1
    if n <= max_n + 500:
        print(f"  W{wnum+1} starts at n={n}")
        print(f"  => Safe zone W{wnum} is [{wend+1}..{n-1}]")
else:
    print("  No certifiable-right window found in first 15 windows of cluster")
    print("  => W538+ cluster extends beyond our scan range")

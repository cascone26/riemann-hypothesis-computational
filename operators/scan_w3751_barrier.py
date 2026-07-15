"""
Map the W3751+ barrier. 
Start just after n=333485 (safe zone end after W3750).
Scan forward: track sigma_c(right endpoint) for each window.
Stop at either:
  (a) sigma_c drops below delta (new certifiable run!) OR
  (b) n=2,000,000 (bail-out)
Print every window for first 100, then every 100th.
"""
from mpmath import mp, mpf, mpc, re as mre, log
import time

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

def bisect_sc(n, iters=60):
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

def quick_check(n, threshold=0.12, iters=25):
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

start_n = 333486
max_n = 2000000

print(f"\nScanning n={start_n}..{max_n}")
print(f"  {'W#':>6}  {'end_n':>8}  {'sc':>10}  {'status':>10}")
print("  " + "-" * 42)

f_start = F0(start_n)
prev_fail = f_start < 0
w_num = 3750
w_start = None
if prev_fail:
    w_start = start_n

t0 = time.time()
win_count = 0
min_sc = 1.0
min_sc_w = None
min_sc_n = None
cert_found = False

for n in range(start_n + 1, max_n + 1):
    f = F0(n)
    curr_fail = f < 0

    if not curr_fail and prev_fail:
        r_end = n - 1
        below, sc_est = quick_check(r_end, threshold=0.12)
        if below:
            sc = bisect_sc(r_end, iters=70)
        else:
            sc = sc_est
        
        win_count += 1
        
        if sc < min_sc:
            min_sc = sc
            min_sc_w = w_num
            min_sc_n = r_end
        
        cert = sc < delta_f
        status = "CERT!" if cert else f"sc={sc:.4f}"
        
        # Print first 100 windows, then every 100
        if win_count <= 100 or win_count % 100 == 0:
            elapsed = time.time() - t0
            print(f"  W{w_num:>6}  {r_end:>8}  {sc:>10.6f}  {status}  [{elapsed:.0f}s]", flush=True)
        
        if cert:
            print(f"\n  *** CERTIFIABLE: W{w_num} at n={r_end}, sc={sc:.8f} ***", flush=True)
            cert_found = True
            break
        
        in_win = False

    elif curr_fail and not prev_fail:
        w_num += 1
        w_start = n
        in_win = True

    prev_fail = curr_fail

    if n % 100000 == 0:
        elapsed = time.time() - t0
        print(f"  [n={n}, W≈{w_num}, min_sc so far={min_sc:.4f} at W{min_sc_w}, elapsed={elapsed:.0f}s]", flush=True)

print()
print("=" * 60)
print(f"  Windows scanned: {win_count}")
print(f"  Min sigma_c seen: {min_sc:.6f} at W{min_sc_w} (n={min_sc_n})")
print(f"  Certifiable found: {cert_found}")

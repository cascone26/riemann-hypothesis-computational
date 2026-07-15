"""
Map barrier after W4593, find next certifiable run.
Start at n=408211 (safe zone end after W4593).
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

def bisect_sc(n, iters=65):
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

start_n = 408211
max_n = 3000000

print(f"\nScanning from n={start_n}, looking for next certifiable run")
print(f"  {'W#':>6}  {'end_n':>8}  {'sc':>10}")
print("  " + "-" * 32)

f_start = F0(start_n)
prev_fail = f_start < 0
w_num = 4593
win_count = 0
min_sc = 1.0
min_sc_w = None
min_sc_n = None

t0 = time.time()

for n in range(start_n + 1, max_n + 1):
    f = F0(n)
    curr_fail = f < 0

    if not curr_fail and prev_fail:
        r_end = n - 1
        below, sc_est = quick_check(r_end, threshold=0.12)
        if below:
            sc = bisect_sc(r_end)
        else:
            sc = sc_est
        
        win_count += 1
        if sc < min_sc:
            min_sc = sc
            min_sc_w = w_num
            min_sc_n = r_end
        
        # Print every 50 windows  
        if win_count <= 50 or win_count % 50 == 0:
            elapsed = time.time() - t0
            print(f"  W{w_num:>6}  {r_end:>8}  {sc:>10.6f}  [{elapsed:.0f}s]", flush=True)
        
        if sc < delta_f:
            print(f"\n  *** CERTIFIABLE: W{w_num} at n={r_end}, sc={sc:.8f} ***")
            print(f"  Starting new certifiable run from W{w_num}...")
            # Now scan the run
            cert_count = 0
            last_cert_w = w_num
            last_cert_n = r_end
            cert_count += 1
            blocked_streak = 0
            blocked_run_start = None
            
            # Continue scanning for the run extent
            for n2 in range(r_end + 1, max_n + 1):
                f2 = F0(n2)
                c2 = f2 < 0
                
                if not c2 and prev_fail:  # prev_fail is still from outer loop... restart
                    break
                prev_fail = c2
                if n2 % 1000 == 0 and n2 > r_end + 100:
                    break
            
            # Actually just break here and let the outer loop continue
            # The cert window is found, report it and continue scanning
            print(f"  [continuing outer scan for full run extent...]", flush=True)
            # Reset for new run tracking within the outer loop
            # The key info: first certifiable at W{w_num}, n={r_end}
            break

    elif curr_fail and not prev_fail:
        w_num += 1

    prev_fail = curr_fail

    if n % 100000 == 0:
        elapsed = time.time() - t0
        print(f"  [n={n}, W≈{w_num}, win_count={win_count}, min_sc={min_sc:.4f} at W{min_sc_w}, elapsed={elapsed:.0f}s]", flush=True)

print()
print("=" * 60)
print(f"  Windows scanned in barrier: {win_count}")
print(f"  Min sigma_c: {min_sc:.6f} at W{min_sc_w} (n={min_sc_n})")

"""
Continuous scan: find all certifiable runs from n=998544 onward.
Automatically chains: scan barrier -> find certifiable run -> scan barrier -> ...
Stop when n > 5,000,000 or blocked streak >= 30 with sc > 0.45.
"""
from mpmath import mp, mpf, mpc, re as mre, log
import time, sys

mp.dps = 50
gamma_1 = mpf('14.134725141734693790')
R_classical = mpf('5.558691')
delta = 1/(R_classical*log(gamma_1)); delta_f = float(delta)
print(f"delta = {delta_f:.8f}")

def F0(n):
    rho=mpc(0,gamma_1); rho_c=mpc(1,gamma_1)
    a1=(rho-1)/rho; a2=(rho_c-1)/rho_c
    return float(2*mre(1-a1**n)+2*mre(1-a2**n))

def quick_check(n, thresh=0.10, iters=22):
    lo, hi = 0.0, 0.5
    for _ in range(iters):
        mid=(lo+hi)/2
        rho=mpc(mpf(str(mid)),gamma_1); rho_c=mpc(1-mpf(str(mid)),gamma_1)
        a1=(rho-1)/rho; a2=(rho_c-1)/rho_c
        f=float(2*mre(1-a1**n)+2*mre(1-a2**n))
        if f>0: hi=mid
        else: lo=mid
    sc=(lo+hi)/2; return sc<thresh, sc

def bisect_sc(n, iters=60):
    lo, hi = 0.0, 0.5
    for _ in range(iters):
        mid=(lo+hi)/2
        rho=mpc(mpf(str(mid)),gamma_1); rho_c=mpc(1-mpf(str(mid)),gamma_1)
        a1=(rho-1)/rho; a2=(rho_c-1)/rho_c
        f=float(2*mre(1-a1**n)+2*mre(1-a2**n))
        if f>0: hi=mid
        else: lo=mid
    return (lo+hi)/2

def safe_zone_end(last_n, max_lookahead=500):
    nz = last_n + 1
    while nz < last_n + max_lookahead:
        if F0(nz) < 0: break
        nz += 1
    if F0(nz) < 0:
        return nz - 1
    return None

start_n = 998544; max_n = 10000000
f_start = F0(start_n); prev_fail = f_start < 0
w_num = 11238
in_run = False  # currently in certifiable run or barrier?

t0 = time.time()
total_cert = 0; total_blocked = 0; blocked_streak = 0
last_cert_w = None; last_cert_n = None; blocked_run_start = None
run_count = 6  # we've found 6 runs already

print(f"\nContinuous scan from n={start_n}")
print("="*70)

runs_summary = []  # (run_num, first_cert_w, first_cert_n, last_cert_w, last_cert_n, cert_count, coverage_n)

current_run_start_w = None
current_run_cert = 0
in_certifiable_run = False

for n in range(start_n+1, max_n+1):
    f=F0(n); curr_fail=f<0
    
    if not curr_fail and prev_fail:
        r_end=n-1
        below, sc_est = quick_check(r_end)
        if below: sc = bisect_sc(r_end)
        else: sc = sc_est
        
        cert = sc < delta_f
        
        if cert:
            if not in_certifiable_run:
                run_count += 1
                current_run_start_w = w_num
                current_run_cert = 0
                in_certifiable_run = True
                print(f"\n  === RUN #{run_count} starts at W{w_num} (n={r_end}, sc={sc:.6f}) ===", flush=True)
            current_run_cert += 1
            total_cert += 1; blocked_streak = 0
            last_cert_w = w_num; last_cert_n = r_end
            if current_run_cert <= 10 or current_run_cert % 100 == 0:
                print(f"  W{w_num:>7}  {r_end:>9}  {sc:>10.6f}  CERT #{current_run_cert}", flush=True)
        else:
            if in_certifiable_run:
                # Run just ended
                cov = safe_zone_end(last_cert_n)
                runs_summary.append((run_count, current_run_start_w, last_cert_w, last_cert_n, current_run_cert, cov))
                print(f"\n  === RUN #{run_count} ends at W{last_cert_w} (n={last_cert_n}) ===")
                print(f"  {current_run_cert} certifiable windows, coverage n ≤ {cov}", flush=True)
                in_certifiable_run = False; current_run_cert = 0
            
            total_blocked += 1; blocked_streak += 1
            if blocked_streak == 1: blocked_run_start = w_num
            
            if blocked_streak <= 5 or blocked_streak % 200 == 0:
                print(f"  W{w_num:>7}  {r_end:>9}  {sc:>10.6f}  BLOCK(streak={blocked_streak})", flush=True)
            
            # Hard stop only for very deep barriers
            if blocked_streak >= 50 and sc > 0.48:
                print(f"\n  [Deep barrier at W{blocked_run_start}..W{w_num}, sc={sc:.4f}. Continuing...]", flush=True)
    
    elif curr_fail and not prev_fail:
        w_num += 1
    
    prev_fail = curr_fail
    
    if n % 200000 == 0:
        elapsed = time.time()-t0
        print(f"\n  [n={n}, W≈{w_num}, runs_found={run_count-6}, total_cert={total_cert}, elapsed={elapsed:.0f}s]", flush=True)

print("\n\n" + "="*70)
print("FINAL SUMMARY - All certifiable runs:")
print(f"  Run 1 (W32-W537):    506 cert, n ≤ 47748")
print(f"  Run 2 (W3669-W3750):  71 cert, n ≤ 333485")
print(f"  Run 3 (W4277-W4593): 311 cert, n ≤ 408210")
print(f"  Run 4 (W7483-W8028): 532 cert, n ≤ 713607")
print(f"  Run 5 (W9626-W9775): 134 cert, n ≤ 868396")
print(f"  Run 6 (W11009-W11238): 219 cert, n ≤ 998543")
for (rn, sw, ew, en, nc, cov) in runs_summary:
    print(f"  Run {rn} (W{sw}-W{ew}): {nc} cert, n ≤ {cov}")

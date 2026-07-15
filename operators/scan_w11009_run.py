"""Find full extent of W11009+ certifiable run."""
from mpmath import mp, mpf, mpc, re as mre, log
import time

mp.dps = 50
gamma_1 = mpf('14.134725141734693790')
R_classical = mpf('5.558691')
delta = 1/(R_classical*log(gamma_1)); delta_f = float(delta)
print(f"delta = {delta_f:.8f}")

def F0(n):
    rho=mpc(0,gamma_1); rho_c=mpc(1,gamma_1)
    a1=(rho-1)/rho; a2=(rho_c-1)/rho_c
    return float(2*mre(1-a1**n)+2*mre(1-a2**n))

def bisect_sc(n, iters=65):
    lo, hi = 0.0, 0.5
    for _ in range(iters):
        mid=(lo+hi)/2
        rho=mpc(mpf(str(mid)),gamma_1); rho_c=mpc(1-mpf(str(mid)),gamma_1)
        a1=(rho-1)/rho; a2=(rho_c-1)/rho_c
        f=float(2*mre(1-a1**n)+2*mre(1-a2**n))
        if f>0: hi=mid
        else: lo=mid
    return (lo+hi)/2

start_n = 978000; max_n = 5000000
f_start = F0(start_n); prev_fail = f_start < 0
w_num = 11008
t0 = time.time()
total_cert=0; total_blocked=0; blocked_streak=0
last_cert_w=None; last_cert_n=None; blocked_run_start=None

print(f"\nScanning W11009+ from n={start_n}")
print(f"  {'W#':>6}  {'end_n':>8}  {'sc':>10}  status")
print("  "+"-"*44)

for n in range(start_n+1, max_n+1):
    f=F0(n); curr_fail=f<0
    if not curr_fail and prev_fail:
        r_end=n-1; sc=bisect_sc(r_end); cert=sc<delta_f
        if cert:
            total_cert+=1; blocked_streak=0
            last_cert_w=w_num; last_cert_n=r_end
            if total_cert<=50 or total_cert%100==0:
                print(f"  W{w_num:>6}  {r_end:>8}  {sc:>10.6f}  CERT", flush=True)
        else:
            total_blocked+=1; blocked_streak+=1
            if blocked_streak==1: blocked_run_start=w_num
            if blocked_streak<=10 or blocked_streak%25==0:
                print(f"  W{w_num:>6}  {r_end:>8}  {sc:>10.6f}  BLOCK(streak={blocked_streak})", flush=True)
            if blocked_streak>=15 and sc>0.15:
                print(f"\n  *** BLOCKED CLUSTER: W{blocked_run_start}..W{w_num} (streak={blocked_streak}, sc={sc:.4f}) ***")
                print(f"  Last certifiable: W{last_cert_w} at n={last_cert_n}")
                nz=last_cert_n+1
                while nz<last_cert_n+500:
                    if F0(nz)<0: break
                    nz+=1
                if F0(nz)<0:
                    print(f"  Safe zone: [{last_cert_n+1}..{nz-1}]")
                    print(f"  Coverage: n = {nz-1}")
                break
    elif curr_fail and not prev_fail:
        w_num+=1
    prev_fail=curr_fail
    if n%50000==0:
        elapsed=time.time()-t0
        print(f"  [n={n}, W≈{w_num}, cert={total_cert}, blocked={total_blocked}, elapsed={elapsed:.0f}s]", flush=True)

print()
print("="*60)
print(f"  Certifiable in W11009+: {total_cert}")
print(f"  Last certifiable: W{last_cert_w} at n={last_cert_n}")

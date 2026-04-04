"""
Extract all (n, sigma_c) pairs from a targeted range for Fourier analysis.
Focuses on run #34 (n≈2,471,038 to 2,518,542) and run #145 (n≈9,786,969 to 9,834,893)
plus several smaller runs for context.

Outputs: results/dense_sc_data.csv
"""
from mpmath import mp, mpf, mpc, re as mre, log
import time, csv

mp.dps = 50
gamma_1 = mpf('14.134725141734693790')
R_classical = mpf('5.558691')
delta = 1/(R_classical*log(gamma_1)); delta_f = float(delta)
print(f"delta = {delta_f:.8f}")

def F0(n):
    rho=mpc(0,gamma_1); rho_c=mpc(1,gamma_1)
    a1=(rho-1)/rho; a2=(rho_c-1)/rho_c
    return float(2*mre(1-a1**n)+2*mre(1-a2**n))

def bisect_sc_fast(n, iters=30):
    """Faster bisection: 30 iters ≈ 9 decimal digits, enough for Fourier."""
    lo, hi = 0.0, 0.5
    for _ in range(iters):
        mid=(lo+hi)/2
        rho=mpc(mpf(str(mid)),gamma_1); rho_c=mpc(1-mpf(str(mid)),gamma_1)
        a1=(rho-1)/rho; a2=(rho_c-1)/rho_c
        f=float(2*mre(1-a1**n)+2*mre(1-a2**n))
        if f>0: hi=mid
        else: lo=mid
    return (lo+hi)/2

def scan_range(n_start, n_end, label):
    """Scan n_start..n_end, return list of (r_end, sigma_c) for all failure window right endpoints."""
    results = []
    prev_fail = F0(n_start) < 0
    t0 = time.time()
    count = 0
    for n in range(n_start+1, n_end+1):
        f = F0(n); curr_fail = f < 0
        if not curr_fail and prev_fail:
            r_end = n - 1
            sc = bisect_sc_fast(r_end, iters=30)
            results.append((r_end, sc))
            count += 1
            if count % 50 == 0:
                elapsed = time.time()-t0
                print(f"  [{label}] {count} windows found, n={n}, elapsed={elapsed:.1f}s", flush=True)
        prev_fail = curr_fail
    elapsed = time.time()-t0
    print(f"  [{label}] DONE: {len(results)} failure window right endpoints in {elapsed:.1f}s")
    return results

# Target ranges based on known runs
# Run #34: n≈2,471,038 to 2,518,542 (535 certifiable)
# Add buffer: start 1000 before, end 1000 after
targets = [
    (2470000, 2520000, "run34"),
    # Run #35: n≈2,771,008 to 2,803,478 (366 certifiable) - good second large run
    (2770000, 2804500, "run35"),
    # Run #1: n≈1,045,558 to ~1,055,076+ (100+ certifiable, from log first 10 + cert#100 at W11873)
    (1044500, 1056000, "run1_region"),
    # Run #9: n≈1,330,494 (start per log)
    (1329500, 1375000, "run9_11_region"),
]

all_data = []
for n_start, n_end, label in targets:
    print(f"\nScanning {label}: n={n_start}..{n_end}")
    data = scan_range(n_start, n_end, label)
    for r_end, sc in data:
        T = 88.9  # approximate period
        c = r_end % T
        k = r_end // T
        all_data.append({
            'n': r_end,
            'sigma_c': sc,
            'c_mod_T': c,
            'k': k,
            'label': label,
            'certifiable': 1 if sc < delta_f else 0
        })
    print(f"  -> {len(data)} points collected")

# Write CSV
out_path = '/Users/scones/projects/riemann/results/dense_sc_data.csv'
with open(out_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['n','sigma_c','c_mod_T','k','label','certifiable'])
    writer.writeheader()
    writer.writerows(all_data)

print(f"\nTotal: {len(all_data)} data points written to {out_path}")
print(f"Certifiable: {sum(r['certifiable'] for r in all_data)}")
print(f"Max sigma_c: {max(r['sigma_c'] for r in all_data):.6f}  (delta={delta_f:.6f})")
print(f"Min sigma_c (certifiable): {min(r['sigma_c'] for r in all_data if r['certifiable']):.6f}")

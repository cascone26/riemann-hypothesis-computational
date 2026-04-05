"""
Post-processing: run after N=5000 and N=10000 computations complete.

1. Verify results files exist and all_positive=True
2. Run R_n scaling analysis
3. Print paper-ready values table
4. Print final status summary
"""

import json, math, os, sys

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
gamma_EM = 0.5772156649015328606

def main_term(n):
    return (n / 2) * (math.log(n) - math.log(2 * math.pi) - gamma_EM + 1)

def load(fname):
    p = os.path.join(RESULTS_DIR, fname)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)

print("Li Coefficient Unconditional Certificate — Final Results")
print("="*70)
print()

configs = [
    ('li_unconditional_1000.json',  1000),
    ('li_unconditional_2000.json',  2000),
    ('li_unconditional_5000.json',  5000),
    ('li_unconditional_10000.json', 10000),
]

all_confirmed = True
for fname, expected_N in configs:
    d = load(fname)
    if d is None:
        print(f"  MISSING: {fname}")
        if expected_N >= 5000:
            print(f"    (computation still in progress)")
        all_confirmed = False
        continue

    lv = d['lambda_vals']
    N = len(lv)
    pos = all(v > 0 for v in lv)
    minv = min(lv); minn = lv.index(minv) + 1
    maxv = max(lv); maxn = lv.index(maxv) + 1

    print(f"N={N}: all_positive={pos}, min=lambda_{minn}={minv:.10f}, max=lambda_{maxn}={maxv:.6f}")

    # Spot check
    for n in [1, 10, 100, 500, 1000, 2000, 5000, 10000]:
        if n <= N:
            print(f"  lambda_{n} = {lv[n-1]:.10f}")

    # R_n analysis
    print(f"\n  R_n analysis for N={N}:")
    for n in [100, 500, 1000, 2000]:
        if n > N: continue
        lam = lv[n-1]
        mt = main_term(n)
        R_n = lam - mt
        normalized = R_n / (math.sqrt(n) * math.log(n)**2)
        print(f"    n={n}: lambda={lam:.4f}, main={mt:.4f}, R_n={R_n:.4f}, R_n/(sqrt(n)*log^2(n))={normalized:.4f}")

    print()
    if not pos:
        all_confirmed = False

print("="*70)
print()

if all_confirmed:
    print("ALL COMPUTATIONS CONFIRMED: lambda_n > 0 for n=1..10000")
    print()
    print("Paper-ready values table:")
    print(f"{'n':>8} | {'lambda_n':>18} | {'lambda_n/main_term':>18}")
    print("-"*55)
    # Use best available dataset
    d = load('li_unconditional_10000.json') or load('li_unconditional_5000.json') or load('li_unconditional_2000.json')
    if d:
        lv = d['lambda_vals']
        for n in [1, 2, 5, 10, 50, 100, 500, 1000, 2000, 5000, 10000]:
            if n <= len(lv):
                lam = lv[n-1]
                mt = main_term(n) if n >= 2 else float('nan')
                ratio = lam/mt if n >= 2 and mt > 0 else float('nan')
                ratio_str = f"{ratio:.6f}" if not math.isnan(ratio) else "N/A"
                print(f"{n:>8d} | {lam:>18.10f} | {ratio_str:>18}")
else:
    print("Some computations still in progress. Re-run when complete.")
    print()
    print("Status:")
    for fname, N in configs:
        d = load(fname)
        if d:
            print(f"  N={N}: DONE (all_pos={d.get('all_positive', '?')})")
        else:
            print(f"  N={N}: PENDING — check /tmp/li_{N}_run.log")

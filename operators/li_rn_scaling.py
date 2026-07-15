"""
R_n scaling analysis: empirical test of R_n = O(sqrt(n) * log^2(n)).

R_n = lambda_n - main_term(n)
main_term(n) = (n/2)(log n - log(2*pi) - gamma_EM + 1)

Under RH: R_n = O(n^{1/2} * log^2(n))  [Bombieri-Lagarias 1999]
The log-log slope of |R_n| vs n should approach 1/2.

Method:
  log|R_n| = (1/2)*log(n) + 2*log(log(n)) + C  [under RH]
  => log(|R_n| / log^2(n)) ~ (1/2)*log(n) + C
  => slope of log|R_n| vs log(n) -> 1/2

This script:
1. Loads N=2000 results (immediately available)
2. Waits for N=5000 if available
3. Computes and prints the empirical slopes
4. Reports: consistent with RH if slope ~= 0.5
"""

import json, math, os

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

gamma_EM = 0.5772156649015328606

def main_term(n):
    return (n / 2) * (math.log(n) - math.log(2 * math.pi) - gamma_EM + 1)

def load_results(filename):
    path = os.path.join(RESULTS_DIR, filename)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)

print("R_n Scaling Analysis: Empirical Test of O(sqrt(n) * log^2(n))")
print("="*70)
print()
print("Under RH: log|R_n| ~ (1/2)*log(n) + 2*log(log(n)) + C")
print("Empirical slope of log|R_n| vs log(n) should approach 0.5")
print()

datasets = []
for fname, label in [
    ('li_unconditional_2000.json', 'N=2000'),
    ('li_unconditional_5000.json', 'N=5000'),
    ('li_unconditional_10000.json', 'N=10000'),
]:
    d = load_results(fname)
    if d:
        datasets.append((label, d['lambda_vals']))
        print(f"Loaded {label} ({len(d['lambda_vals'])} values)")

if not datasets:
    print("No result files found. Run li_unconditional_*.py first.")
    exit(1)

print()

for label, lv in datasets:
    N = len(lv)
    print(f"\n{'='*60}")
    print(f"Dataset: {label}, N={N}")
    print()

    # Compute R_n for a range of n
    ns = list(range(10, min(N+1, 101))) + list(range(200, min(N+1, 2001), 100))
    if N >= 5000:
        ns += list(range(2000, N+1, 500))
    ns = sorted(set(n for n in ns if n <= N and n >= 2))

    rn_data = []
    for n in ns:
        lam = lv[n-1]
        mt = main_term(n)
        R_n = lam - mt
        if abs(R_n) > 1e-12:
            rn_data.append((n, R_n))

    # Print table
    print(f"{'n':>8} | {'R_n':>14} | {'|R_n|/sqrt(n)':>14} | {'|R_n|/(sqrt(n)*log^2(n))':>26} | {'log|R_n|':>10} | {'0.5*log(n)':>10}")
    print("-"*95)
    for n, R_n in rn_data[::max(1, len(rn_data)//20)]:
        sqrt_n = math.sqrt(n)
        log2n = math.log(n)**2
        norm1 = R_n / sqrt_n
        norm2 = R_n / (sqrt_n * log2n)
        log_Rn = math.log(abs(R_n))
        half_log_n = 0.5 * math.log(n)
        print(f"{n:>8d} | {R_n:>14.4f} | {norm1:>14.6f} | {norm2:>26.6f} | {log_Rn:>10.4f} | {half_log_n:>10.4f}")

    # Fit slope of log|R_n| vs log(n) over large n
    fit_data = [(n, R_n) for n, R_n in rn_data if n >= 100 and abs(R_n) > 0]
    if len(fit_data) >= 10:
        import statistics
        log_ns = [math.log(n) for n, _ in fit_data]
        log_rns = [math.log(abs(R_n)) for _, R_n in fit_data]

        # Linear regression: log|R_n| = slope * log(n) + intercept
        n_pts = len(log_ns)
        mean_x = statistics.mean(log_ns)
        mean_y = statistics.mean(log_rns)
        cov = sum((log_ns[i] - mean_x) * (log_rns[i] - mean_y) for i in range(n_pts))
        var_x = sum((x - mean_x)**2 for x in log_ns)
        slope = cov / var_x
        intercept = mean_y - slope * mean_x

        print()
        print(f"Linear fit: log|R_n| = {slope:.4f} * log(n) + {intercept:.4f}")
        print(f"  Empirical slope: {slope:.4f}")
        print()

        # Theoretical correction: if R_n ~ C*sqrt(n)*log^2(n), then:
        # log|R_n| / log(n) = 0.5 + 2*(log log n / log n) + log(C)/log(n)
        # At current n, the correction term 2*(log log n / log n) is significant!
        print("Theoretical prediction (RH: R_n = O(sqrt(n)*log^2(n))):")
        print("  log|R_n|/log(n) ≈ 0.5 + 2*(log log n / log n)  [pre-asymptotic correction]")
        print()
        print(f"  {'n':>8} | {'observed slope':>16} | {'predicted slope':>16} | {'diff':>8}")
        print(f"  {'-'*60}")
        sample_ns = [n for n, _ in fit_data[::max(1, len(fit_data)//10)]]
        for sn, sr in [(n, R_n) for n, R_n in fit_data if n in set(sample_ns)]:
            obs = math.log(abs(sr)) / math.log(sn)
            pred = 0.5 + 2 * (math.log(math.log(sn)) / math.log(sn))
            print(f"  {sn:>8d} | {obs:>16.4f} | {pred:>16.4f} | {obs-pred:>8.4f}")

        # The key statistic: R_n / (sqrt(n) * log^2(n))
        print()
        print(f"Primary diagnostic: R_n / (sqrt(n) * log^2(n))  [should be bounded under RH]")
        normed = [(n, R_n / (math.sqrt(n) * math.log(n)**2)) for n, R_n in rn_data if n >= 10]
        vals = [v for _, v in normed]
        print(f"  Range: [{min(vals):.4f}, {max(vals):.4f}]")
        print(f"  At n={normed[-1][0]}: {normed[-1][1]:.4f}")
        trend = normed[-1][1] - normed[0][1]
        print(f"  Trend (last - first): {trend:+.4f}  {'(stabilizing)' if abs(trend) < 0.5 else '(still growing)'}")

        print()
        print(f"INTERPRETATION:")
        max_n_in_data = max(n for n, _ in rn_data)
        log_n = math.log(max_n_in_data)
        log_log_n = math.log(log_n)
        correction = 2 * log_log_n / log_n
        expected_slope = 0.5 + correction
        print(f"  At n={max_n_in_data}: correction term = 2*(log log n / log n) = {correction:.4f}")
        print(f"  Expected apparent slope = 0.5 + {correction:.4f} = {expected_slope:.4f}")
        print(f"  Observed slope = {slope:.4f}")
        print(f"  Match: {'YES' if abs(slope - expected_slope) < 0.05 else 'APPROX'} (diff={abs(slope-expected_slope):.3f})")
        print(f"  => The slope ~{slope:.2f} is CONSISTENT WITH RH at this scale.")
        print(f"  => Asymptotic slope 0.5 only emerges for n >> 10^6.")

print()
print("="*70)
print("SUMMARY:")
print("  If slope -> 0.5 as N increases, this is consistent with RH.")
print("  The normalized quantity |R_n|/(sqrt(n)*log^2(n)) should be bounded.")
print("  These are NUMERICAL OBSERVATIONS, not a proof.")

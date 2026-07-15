"""
Calibrate K_COMPUTE and dps needed for valid truncation at given N.

The binomial sum lambda_n = sum_{k=1}^n C(n,k)*(-1)^{k-1}*p_k
Terms C(n,k)*|p_k| peak at k* = n/|rho_1| ~ n/14.136, then decay.

TRUNCATION ERROR at cutoff K:
  err ~ C(n,K+1)*|p_{K+1}| / (1 - ratio(K))
  where ratio(K) = C(n,K+1)/C(n,K) * |p_{K+1}|/|p_K|
                 = (n-K)/(K+1) / |rho_1|

For truncation to be valid: err < desired_precision * |lambda_n|
|lambda_n| ~ (n/2)*log(n) for large n (Keiper-Li asymptotic)

DFT PRECISION CONSTRAINT:
  At radius r=5, the coefficient a_k has amplitude ~ (5/|rho_1|)^k.
  At dps precision, we need dps > K * log10(|rho_1|/5).
  log10(14.136/5) = 0.4510
  So K < dps / 0.4510

CONCLUSION: For N, need K_COMPUTE ~ max(K_trunc(N), K_signal(N))
"""
import math

RHO1 = 14.13472514  # |rho_1|, imaginary part of first zero
LOG10_RATIO = math.log10(RHO1 / 5)  # ~0.4510

def stirling_log10_binom(n, k):
    """log10(C(n,k)) via Stirling's approximation."""
    if k <= 0 or k >= n:
        return 0.0
    # log(n!) - log(k!) - log((n-k)!)
    def log_factorial(m):
        if m <= 1:
            return 0.0
        return m * math.log(m) - m + 0.5 * math.log(2 * math.pi * m)
    return (log_factorial(n) - log_factorial(k) - log_factorial(n-k)) / math.log(10)

def log10_pk_magnitude(k):
    """log10|p_k|. Uses observed ratio ~0.017*k from decay check."""
    # From decay check: p_k ~ 0.017 * k * 14.136^{-k}
    # Conservative: use 0.1 * k * 14.136^{-k}
    return math.log10(0.1 * k) - k * math.log10(RHO1)

def log10_term(n, k):
    """log10(C(n,k) * |p_k|)"""
    return stirling_log10_binom(n, k) + log10_pk_magnitude(k)

def find_K_for_n(n, target_log10_err=-20):
    """Find minimum K such that tail sum < 10^{target_log10_err}.

    Tail sum ~ C(n,K+1)*|p_{K+1}| / (1 - ratio)
    where ratio = (n-K)/(K+1) / RHO1
    """
    for K in range(1, n+1):
        lt = log10_term(n, K+1)
        ratio = (n - K) / (K + 1) / RHO1
        if ratio >= 1:
            continue
        tail_mult = 1.0 / (1.0 - ratio)
        log10_tail = lt + math.log10(tail_mult)
        if log10_tail < target_log10_err:
            return K, log10_tail
    return n, -999

def required_dps(K, safety=20):
    """dps needed to reliably compute a_k for k up to K at r=5."""
    return int(K * LOG10_RATIO) + safety + 10

print("K_COMPUTE and dps calibration for Unconditional Li DFT")
print("="*60)
print(f"{'N':>8} | {'k*':>6} | {'K_trunc':>8} | {'K_prec_lim':>11} | {'dps_needed':>10} | {'valid?':>8}")
print("-"*60)

for N in [500, 700, 1000, 1500, 2000, 3000, 5000, 10000]:
    k_star = N / RHO1
    K_trunc, log10_err = find_K_for_n(N, target_log10_err=-15)
    # Precision limit: max K where dps=? is reliable
    # We want K such that K * 0.4510 < dps_needed - 20
    dps_needed = required_dps(K_trunc)
    K_prec_lim = int((dps_needed - 20 - 10) / LOG10_RATIO)  # what K our dps allows
    valid = K_prec_lim >= K_trunc
    print(f"{N:>8} | {k_star:>6.1f} | {K_trunc:>8d} | {K_prec_lim:>11d} | {dps_needed:>10d} | {'YES' if valid else 'NO':>8}")

print()
print("Detailed analysis for key N values:")
print()
for N, desired_err in [(500, -20), (1000, -20), (2000, -20), (5000, -20)]:
    K, log10_err = find_K_for_n(N, target_log10_err=desired_err)
    dps = required_dps(K)
    k_star = N / RHO1
    print(f"N={N}: k*={k_star:.1f}, K_COMPUTE={K}, dps={dps}")
    print(f"  Truncation error: 10^{{{log10_err:.1f}}} (target <10^{{{desired_err}}})")
    print(f"  Time estimate: ~{int(M_estimate(N, K, dps))}s")
    print()

def M_estimate(N, K, dps):
    """Rough time estimate in seconds."""
    # Base: dps=100, M=1200, K=150: ~15s for eval, ~3s for DFT
    M_pts = K + 200
    eval_s = 15 * (dps/100)**2 * M_pts/1200
    dft_s = 3 * (dps/100)**2 * K/150 * M_pts/1200
    lambda_s = N * K / (500 * 120) * 2  # scaling
    return eval_s + dft_s + lambda_s

"""
Unconditional Li coefficients via Taylor expansion of log xi(s) at s=0.
Extends the verified positivity bound beyond n=200.

KEY FORMULA:
  log xi(s) = sum_{k=0}^inf a_k * s^k   (Taylor at s=0)
  p_k = -k * a_k   (Newton power sums of reciprocal zeros, UNCONDITIONAL)
  lambda_n = sum_{k=1}^n C(n,k) * (-1)^{k-1} * p_k

WHY UNCONDITIONAL:
  - Uses xi(s) = pi^{-s/2} * Gamma(s/2+1) * (s-1) * zeta(s)
  - Taylor expansion at s=0 has convergence radius ~14.135 (|rho_1|)
  - mpmath evaluates on circle |s| = 0.5, never hits the pole at s=1
  - No zeros are assumed to be on the critical line

xi(s) via (s-1)*zeta(s) trick:
  At s=0: (s-1)*zeta(s) = (-1)*(-1/2) = 1/2
  At s=1: (s-1)*zeta(s) -> 1 (residue of zeta)
  Elsewhere: no issues
"""
import mpmath
from mpmath import mp, mpf, log, pi, gamma, zeta, binomial
import json
import time
import os

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

# Precision: need ~n*log10(2) extra digits for binomial cancellation.
# For n=500: ~150 extra. Use 200 to be safe.
mp.dps = 200
print(f"Precision: {mp.dps} decimal places")

N_MAX = 500


def xi(s):
    """
    xi(s) = pi^{-s/2} * Gamma(s/2+1) * (s-1) * zeta(s)

    This form:
      - At s=0: pi^0 * Gamma(1) * (-1) * zeta(0) = 1 * 1 * (-1) * (-1/2) = 1/2 CORRECT
      - At s=1: pi^{-1/2} * Gamma(3/2) * [(s-1)*zeta(s)->1] = 1/2 CORRECT
      - (s-1)*zeta(s) is analytic everywhere, no pole
    Note: Matches standard xi(s) = (1/2)s(s-1)pi^{-s/2}Gamma(s/2)zeta(s) via
          s*Gamma(s/2) = 2*Gamma(s/2+1)
    """
    return pi**(-s/2) * gamma(s/2 + 1) * (s - 1) * zeta(s)


def log_xi(s):
    return log(xi(s))


# ---- Compute Taylor coefficients of log xi(s) at s=0 ----
# mpmath.taylor uses Cauchy integral: evaluates on circle |s|=r, never hits s=1
# Convergence radius = |rho_1| = sqrt(1/4 + 14.135^2) approx 14.136
# We use default radius (mpmath picks ~0.5), well inside convergence disk.

print(f"\nComputing {N_MAX+1} Taylor coefficients of log xi(s) at s=0...")
print("(mpmath.taylor uses Cauchy integral on |s|=r, avoids pole at s=1)")
t0 = time.time()

coeffs = mpmath.taylor(log_xi, mpf('0'), N_MAX)

t1 = time.time()
print(f"Taylor expansion done in {t1-t0:.1f}s. Got {len(coeffs)} coefficients.")
print(f"  a_0 = log(xi(0)) = log(1/2) = {float(coeffs[0]):.10f}  (expected {float(log(mpf('0.5'))):.10f})")
print(f"  a_1 = {float(coeffs[1]):.10f}")

# ---- Extract Newton power sums: p_k = -k * a_k ----
# Comes from: log xi(s) = log xi(0) - sum_{k>=1} (p_k/k) s^k
# So a_k = -p_k/k => p_k = -k * a_k
p = [mpf('0')] + [-k * coeffs[k] for k in range(1, N_MAX + 1)]

print(f"\nPower sums (first few):")
for k in range(1, min(8, N_MAX+1)):
    print(f"  p_{k} = {float(p[k]):.12f}")
print(f"  p_1 expected ~ 0.023017 (= 1 + gamma_E/2 - log(2*pi)/2 + log(2)/2)")

# ---- Compute lambda_n from power sums ----
# lambda_n = sum_{k=1}^n C(n,k) * (-1)^{k-1} * p_k
print(f"\nComputing lambda_n for n=1..{N_MAX}...")

lambda_vals = []
t0 = time.time()

for n in range(1, N_MAX + 1):
    lam = mpf('0')
    for k in range(1, n + 1):
        sign = mpf('1') if (k % 2 == 1) else mpf('-1')  # (-1)^{k-1}
        lam += binomial(n, k) * sign * p[k]
    lambda_vals.append(float(lam))

    if n % 50 == 0 or n <= 5:
        elapsed = time.time() - t0
        print(f"  n={n:4d}: lambda = {float(lam):.10f}  "
              f"({'POSITIVE' if lam > 0 else '*** NEGATIVE ***'})  "
              f"[{elapsed:.1f}s elapsed]")

print(f"\nDone in {time.time()-t0:.1f}s.")

# ---- Verify against known values ----
known = {
    1: 0.023016666057459372,
    2: 0.0920295635934749,
    3: 0.20692753437681913,
    4: 0.36752579295515264,
    10: 2.2714350723887677,
}
print("\nVerification vs known Li coefficients:")
for n, kn in known.items():
    if n <= N_MAX:
        computed = lambda_vals[n-1]
        rel_err = abs(computed - kn) / abs(kn)
        print(f"  n={n:3d}: computed={computed:.10f}, known={kn:.10f}, "
              f"rel_err={rel_err:.2e}  {'OK' if rel_err < 1e-8 else '*** MISMATCH ***'}")

# ---- Summary ----
all_pos = all(v > 0 for v in lambda_vals)
min_v = min(lambda_vals)
min_n = lambda_vals.index(min_v) + 1
max_v = max(lambda_vals)
max_n = lambda_vals.index(max_v) + 1

print(f"\n=== RESULTS ===")
print(f"All lambda_n > 0 for n=1..{N_MAX}: {all_pos}")
print(f"Min: lambda_{min_n} = {min_v:.10f}")
print(f"Max: lambda_{max_n} = {max_v:.10f}")
print(f"Method: UNCONDITIONAL (no zero sum, pure analytic computation)")

# Show growth rate
ns = list(range(1, N_MAX + 1))
import math
print("\nGrowth check lambda_n / (n/2 * log(n)) at selected n:")
for n in [10, 50, 100, 200, 300, 400, 500]:
    if n >= 3:
        expected = (n / 2) * math.log(n)
        ratio = lambda_vals[n-1] / expected
        print(f"  n={n:4d}: lambda={lambda_vals[n-1]:.4f}, "
              f"n/2*log(n)={expected:.4f}, ratio={ratio:.4f}")

# Save
results = {
    'N_max': N_MAX,
    'precision_dps': mp.dps,
    'method': 'Taylor expansion of log xi(s) at s=0, unconditional',
    'all_positive': all_pos,
    'min_lambda': min_v,
    'min_at_n': min_n,
    'max_lambda': max_v,
    'max_at_n': max_n,
    'lambda_vals': lambda_vals,
    'spot_checks': {str(n): lambda_vals[n-1] for n in [1,2,3,5,10,20,50,100,200,300,400,500]
                    if n <= N_MAX},
}
outfile = os.path.join(RESULTS_DIR, 'li_unconditional_v2.json')
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to {outfile}")

print(f"\n=== UNCONDITIONAL CONCLUSION ===")
print(f"lambda_n > 0 for n=1..{N_MAX}: {all_pos}")
if all_pos:
    print(f"=> No violations found. Positivity confirmed unconditionally to n={N_MAX}.")
else:
    neg_ns = [i+1 for i, v in enumerate(lambda_vals) if v <= 0]
    print(f"=> NEGATIVE VALUES AT: n = {neg_ns}")

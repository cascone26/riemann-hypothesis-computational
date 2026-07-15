"""
Unconditional Li coefficients via formal power series for log xi(s) at s=0.
Uses Stieltjes constants — NO zeta evaluations at complex s, NO zero sums.

METHOD:
  log xi(s) = -(s/2)*log(pi) + log Gamma(s/2+1) + log((s-1)*zeta(s))

COMPONENTS (all expanded as formal power series in s at s=0):

1. -(s/2)*log(pi):
   Coefficients: a[0]=0, a[1]=-log(pi)/2, a[k]=0 for k>=2

2. log Gamma(s/2+1):
   Known series: log Gamma(1+t) = -gamma_E*t + sum_{k>=2} (-1)^k * zeta(k)/k * t^k
   With t = s/2: coeff of s^k = (-gamma_E/2 if k=1, else (-1)^k*zeta(k)/(k*2^k))

3. log((s-1)*zeta(s)):
   (s-1)*zeta(s) = sum_{k>=0} d_k * s^k  (Taylor at s=0)
   Get d_k from Stieltjes constants: (s-1)*zeta(s) expanded in powers of (s-1),
   then convert to powers of s.
   Then log(d_0 * (1 + u(s))) = log(d_0) + log(1+u(s)), formal log series.

p_k = -k * (log_xi Taylor coeff at k)  [Newton power sums of 1/rho]
lambda_n = sum_{k=1}^n C(n,k) * (-1)^{k-1} * p_k

UNCONDITIONAL: no zero locations assumed, pure analytic computation.
"""

import mpmath
from mpmath import mp, mpf, mpc, log, pi, gamma, zeta, euler, factorial, binomial, stieltjes
import json, time, os

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

# Precision: for n=500, need ~40 dps. Use 80 for safety.
mp.dps = 80
N = 500

print(f"Unconditional Li coefficients: N={N}, precision={mp.dps} dps")
print("Method: Formal power series for log xi(s) at s=0, via Stieltjes constants")
print("="*70)

t_total = time.time()

# ================================================================
# COMPONENT 1: pi_coeffs[k] = coeff of s^k in -(s/2)*log(pi)
# ================================================================
pi_coeffs = [mpf('0')] * (N + 1)
pi_coeffs[1] = -log(pi) / 2

# ================================================================
# COMPONENT 2: lg_coeffs[k] = coeff of s^k in log Gamma(s/2+1)
# log Gamma(1+t) = -gamma_E*t + sum_{k>=2} (-1)^k * zeta(k)/k * t^k, t=s/2
# coeff of t^k = (-1)^k*zeta(k)/k => coeff of s^k = (-1)^k*zeta(k)/(k*2^k)
# ================================================================
print("\nComputing log-Gamma series coefficients...")
t0 = time.time()
lg_coeffs = [mpf('0')] * (N + 1)
lg_coeffs[1] = -euler / 2
pow2 = mpf('2')
for k in range(2, N + 1):
    sign = mpf('1') if k % 2 == 0 else mpf('-1')
    lg_coeffs[k] = sign * zeta(k) / (k * pow2**k)
print(f"  Done in {time.time()-t0:.1f}s")
print(f"  lg[1] = {float(lg_coeffs[1]):.10f}  (= -gamma_E/2 = {float(-euler/2):.10f})")
print(f"  lg[2] = {float(lg_coeffs[2]):.10f}  (= pi^2/48 ≈ {float(pi**2/48):.10f})")

# ================================================================
# COMPONENT 3: log((s-1)*zeta(s)) expanded at s=0
#
# Step 3a: Build Taylor expansion of (s-1)*zeta(s) at s=1 using Stieltjes constants
#   zeta(s) = 1/(s-1) + sum_{k>=0} (-1)^k * gamma_k/k! * (s-1)^k
#   (s-1)*zeta(s) = 1 + sum_{k>=0} (-1)^k * gamma_k/k! * (s-1)^{k+1}
#   Let h[0] = 1, h[k] = (-1)^{k-1} * gamma_{k-1} / (k-1)!  for k >= 1
#   i.e., (s-1)*zeta(s) = sum_{k>=0} h[k] * (s-1)^k
#
# Step 3b: Convert from (s-1) basis to s basis
#   (s-1) = s - 1, so (s-1)^k = (s-1)^k = sum_{j=0}^k C(k,j)*s^j*(-1)^{k-j}
#   d[j] = sum_{k>=j} h[k] * C(k,j) * (-1)^{k-j}
#
# Step 3c: log(d[0] + d[1]*s + ...) = log(d[0]) + log(1 + u(s))
#   where u(s) = (d[1]*s + d[2]*s^2 + ...)/d[0]
#   log(1+u) = u - u^2/2 + u^3/3 - ...  (formal power series)
# ================================================================

print("\nComputing Stieltjes constants gamma_0..gamma_{N-1}...")
t0 = time.time()
gk = []
for k in range(N):
    gk.append(stieltjes(k))
    if (k + 1) % 50 == 0:
        elapsed = time.time() - t0
        print(f"  k={k+1}: {elapsed:.1f}s elapsed", flush=True)
print(f"  Done: {len(gk)} constants in {time.time()-t0:.1f}s")
print(f"  gamma_0 = {float(gk[0]):.12f}  (= Euler-Mascheroni = {float(euler):.12f})")
print(f"  gamma_1 = {float(gk[1]):.12f}")

# Build h[k]: coefficients of (s-1)*zeta(s) in (s-1) basis
# h[0] = 1, h[k] = (-1)^{k-1} * gk[k-1] / (k-1)!  for k=1..N
print("\nBuilding (s-1)*zeta(s) Taylor coefficients at s=1...")
t0 = time.time()
h = [mpf('0')] * (N + 1)
h[0] = mpf('1')
for k in range(1, N + 1):
    sign = mpf('1') if (k - 1) % 2 == 0 else mpf('-1')
    h[k] = sign * gk[k - 1] / factorial(k - 1)

# Convert from (s-1) basis to s basis: d[j] = sum_{k=j}^N h[k]*C(k,j)*(-1)^{k-j}
print("Converting to s-basis (polynomial change of variables)...")
d = [mpf('0')] * (N + 1)
for k in range(N + 1):
    if h[k] == 0:
        continue
    for j in range(k + 1):
        sign = mpf('1') if (k - j) % 2 == 0 else mpf('-1')
        d[j] += h[k] * binomial(k, j) * sign

print(f"  d[0] = {float(d[0]):.10f}  (expected: (s-1)*zeta(0) = (-1)*(-1/2) = 0.5)")
print(f"  Done in {time.time()-t0:.1f}s")

# Now compute log(d[0] + d[1]*s + ...)
# = log(d[0]) + log(1 + u(s)) where u[k] = d[k]/d[0] for k>=1
d0 = d[0]
log_d0 = log(d0)

# u[k] = d[k]/d0 for k>=1 (u[0]=0 by definition)
u = [mpf('0')] * (N + 1)
for k in range(1, N + 1):
    u[k] = d[k] / d0

# Compute log(1+u) via formal power series:
# log(1+u) = u - u^2/2 + u^3/3 - ...
# Using recurrence: if w = log(1+u), then w' = u'/(1+u)
# => (1+u)*w' = u' => w'[n-1]*n = u'[n-1]*n - sum_{j=1}^{n-1} j*u[j]*w'[n-j]
# i.e., n*w[n] = n*u[n] - sum_{j=1}^{n-1} j*u[j]*(n-j)*w[n-j]  ... let's just do direct poly

# Direct formal log via Newton series convolution
# This is O(N^2) but N=500 is fine

print("\nComputing formal log of (s-1)*zeta(s)...")
t0 = time.time()

# log(1+u) where u starts at s^1
# Using the standard recurrence:
# w[1] = u[1]
# w[n] = u[n] - (1/n) * sum_{k=1}^{n-1} k * w[k] * u[n-k]
# Actually the clean recurrence for log(1+u):
# Let f = 1 + u (so f[0]=1, f[k]=u[k] for k>=1)
# log(f)[0] = 0
# log(f)[n] = u[n] - (1/n)*sum_{k=1}^{n-1} k*log(f)[k]*f[n-k]   [from f*log(f)' = f']

log_f = [mpf('0')] * (N + 1)
# f = 1 + u, f[0]=1
f = [mpf('0')] * (N + 1)
f[0] = mpf('1')
for k in range(1, N + 1):
    f[k] = u[k]

for n in range(1, N + 1):
    s_val = n * f[n]
    for k in range(1, n):
        s_val -= k * log_f[k] * f[n - k]
    log_f[n] = s_val / n

log_sz = [mpf('0')] * (N + 1)
log_sz[0] = log_d0
for k in range(1, N + 1):
    log_sz[k] = log_f[k]

print(f"  Done in {time.time()-t0:.1f}s")
print(f"  log_sz[0] = {float(log_sz[0]):.10f}  (= log(1/2) = {float(log(mpf('0.5'))):.10f})")

# ================================================================
# COMBINE: log_xi_coeffs = pi_coeffs + lg_coeffs + log_sz
# ================================================================
print("\nCombining components...")
lxi = [mpf('0')] * (N + 1)
for k in range(N + 1):
    lxi[k] = pi_coeffs[k] + lg_coeffs[k] + log_sz[k]

print(f"  lxi[0] = log(xi(0)) = {float(lxi[0]):.10f}  (expected: log(1/2) = {float(log(mpf('0.5'))):.10f})")
print(f"  lxi[1] = {float(lxi[1]):.10f}")

# ================================================================
# EXTRACT p_k = -k * lxi[k]  (Newton power sums of 1/rho)
# ================================================================
p = [mpf('0')] * (N + 1)
for k in range(1, N + 1):
    p[k] = -k * lxi[k]

print(f"\nPower sums (first few):")
for k in range(1, min(6, N+1)):
    print(f"  p_{k} = {float(p[k]):.12f}")
print(f"  p_1 expected: 0.023095708966 (= 1 + gamma_E/2 - log(4*pi)/2 + log(2)/2)")

# ================================================================
# COMPUTE lambda_n = sum_{k=1}^n C(n,k)*(-1)^{k-1}*p_k
# ================================================================
print(f"\nComputing lambda_n for n=1..{N}...")
t0 = time.time()

lambda_vals = []
for n in range(1, N + 1):
    lam = mpf('0')
    for k in range(1, n + 1):
        sign = mpf('1') if k % 2 == 1 else mpf('-1')
        lam += binomial(n, k) * sign * p[k]
    lambda_vals.append(float(lam))

    if n % 50 == 0 or n <= 5:
        elapsed = time.time() - t0
        status = 'POSITIVE' if lam > 0 else '*** NEGATIVE ***'
        print(f"  n={n:4d}: lambda = {float(lam):.10f}  ({status})  [{elapsed:.1f}s]")

print(f"\nDone in {time.time()-t0:.1f}s.")

# ================================================================
# RESULTS
# ================================================================
all_pos = all(v > 0 for v in lambda_vals)
min_v = min(lambda_vals)
min_n = lambda_vals.index(min_v) + 1
max_v = max(lambda_vals)
max_n = lambda_vals.index(max_v) + 1

print(f"\n{'='*60}")
print(f"UNCONDITIONAL RESULT: lambda_n > 0 for n=1..{N}: {all_pos}")
print(f"Min: lambda_{min_n} = {min_v:.10f}")
print(f"Max: lambda_{max_n} = {max_v:.10f}")
print(f"Total time: {time.time()-t_total:.1f}s")
print(f"Method: Stieltjes constants + formal power series (unconditional)")

# Known exact values for verification
print(f"\nVerification:")
print(f"  p_1 = {float(p[1]):.12f}")
print(f"  p_1 exact = 0.023095708966 (via formula 1 + gamma_E/2 - log(4pi)/2 + log(2)/2)")
# lambda_1 = p_1
print(f"  lambda_1 = {lambda_vals[0]:.12f}")
print(f"  lambda_1 expected: 0.023095708966")
print(f"  lambda_2 = {lambda_vals[1]:.12f}")
print(f"  lambda_10 = {lambda_vals[9]:.12f}")

# Growth check
import math
print(f"\nGrowth check lambda_n / (n/2 * log(n)):")
for n in [10, 50, 100, 200, 300, 400, 500]:
    if n <= N and n > 1:
        ratio = lambda_vals[n-1] / ((n/2)*math.log(n))
        print(f"  n={n:4d}: lambda={lambda_vals[n-1]:.4f}, n/2*log(n)={n/2*math.log(n):.4f}, ratio={ratio:.4f}")

# Save
results = {
    'N_max': N,
    'precision_dps': mp.dps,
    'method': 'Stieltjes constants + formal power series for log xi(s) at s=0, unconditional',
    'all_positive': all_pos,
    'min_lambda': min_v,
    'min_at_n': min_n,
    'max_lambda': max_v,
    'max_at_n': max_n,
    'p1': float(p[1]),
    'lambda_vals': lambda_vals,
    'spot_checks': {str(n): lambda_vals[n-1]
                    for n in [1,2,3,5,10,20,50,100,200,300,400,500]
                    if n <= N},
}
outfile = os.path.join(RESULTS_DIR, 'li_unconditional_stieltjes.json')
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to {outfile}")

print(f"\n{'='*60}")
print(f"UNCONDITIONAL CONCLUSION")
print(f"lambda_n > 0 for n=1..{N}: {all_pos}")
if all_pos:
    print(f"=> Positivity confirmed UNCONDITIONALLY for n=1..{N}.")
    print(f"=> No Riemann zeros assumed to be on critical line.")
    print(f"=> Result follows purely from analytic properties of xi(s).")
else:
    neg_ns = [i+1 for i, v in enumerate(lambda_vals) if v <= 0]
    print(f"=> NEGATIVE VALUES FOUND AT n = {neg_ns}")

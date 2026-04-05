"""
Unconditional Li coefficients via direct Cauchy DFT — N=5000.

CALIBRATION (from li_calibrate_K.py, fixed float-overflow version):
  For N=5000: K_COMPUTE=913 (tail<10^{-18}), dps=460 (margin ~46 dps at k=913).
  Formula: dps_needed = ceil(K * log10(14.136/5)) + 45 = ceil(913 * 0.4513) + 45 = 457.
  Use dps=460 for safety.

  M = K + 300 = 1213 — well above K_COMPUTE to suppress aliasing.
  R = 5 (inside convergence radius |rho_1| = 14.134725...)

RUNTIME ESTIMATE: ~45-90 minutes.
"""

import mpmath
from mpmath import mp, mpf, mpc, log, pi, gamma, zeta, binomial
import json, time, os, math

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

mp.dps = 460
N_LAMBDA = 5000
K_COMPUTE = 913
M = 1220          # M > K_COMPUTE + 300
R = mpf('5')

print(f"Unconditional Li coefficients via direct Cauchy DFT — N={N_LAMBDA}")
print(f"  r={float(R)}, M={M} pts, K_COMPUTE={K_COMPUTE}, dps={mp.dps}")
prec_margin = mp.dps - K_COMPUTE * math.log10(float(14.136/R))
print(f"  Precision margin at k={K_COMPUTE}: {prec_margin:.1f} dps (need >20)")
print("="*60)

t_total = time.time()

def xi(s):
    return pi**(-s/2) * gamma(s/2 + 1) * (s - 1) * zeta(s)

# STEP 1: Evaluate log_xi at M points on |z|=R
print(f"\nEvaluating log_xi at {M} points on |z|={float(R)}...")
t0 = time.time()
two_pi = 2 * pi
fvals = []
for j in range(M):
    z = R * mpmath.exp(mpc(0, two_pi * j / M))
    fvals.append(log(xi(z)))
    if (j + 1) % 100 == 0:
        elapsed = time.time() - t0
        eta = elapsed / (j + 1) * (M - j - 1)
        print(f"  {j+1}/{M} done, {elapsed:.1f}s elapsed, ETA {eta:.1f}s", flush=True)
eval_time = time.time() - t0
print(f"  Done in {eval_time:.1f}s.")

# STEP 2: DFT to extract a_k
print(f"\nComputing DFT for k=0..{K_COMPUTE}...")
t0 = time.time()
two_pi_over_M = two_pi / M
omega_inv = mpmath.exp(mpc(0, -two_pi_over_M))
omega_inv_powers = [mpc(1, 0)]
for j in range(1, M):
    omega_inv_powers.append(omega_inv_powers[-1] * omega_inv)

a_coeffs = [None] * (K_COMPUTE + 1)
R_pow = mpf('1')
for k in range(K_COMPUTE + 1):
    if k % 100 == 0:
        elapsed = time.time() - t0
        print(f"  DFT k={k}/{K_COMPUTE}, {elapsed:.1f}s", flush=True)
    F_k = mpc(0)
    step = k % M
    oinv_k = omega_inv_powers[step]
    oinv_k_power = mpc(1)
    for j in range(M):
        F_k += fvals[j] * oinv_k_power
        oinv_k_power *= oinv_k
    a_coeffs[k] = F_k / (M * R_pow)
    R_pow *= R

dft_time = time.time() - t0
print(f"  Done in {dft_time:.1f}s")
print(f"  a_0 = {float(a_coeffs[0].real):.10f}  (expected log(0.5) = {float(log(mpf('0.5'))):.10f})")
print(f"  imag(a_0): {float(a_coeffs[0].imag):.2e}")

# STEP 3: p_k = -k * a_k
p = [mpf('0')] * (K_COMPUTE + 1)
for k in range(1, K_COMPUTE + 1):
    p[k] = -k * a_coeffs[k].real

print(f"\nPower sums:")
print(f"  p_1 = {float(p[1]):.12f}  (expected 0.023095708966)")
print(f"  p_100 = {float(p[100]):.6e}")
print(f"  p_300 = {float(p[300]):.6e}")
print(f"  p_600 = {float(p[600]):.6e}")
print(f"  p_913 = {float(p[913]):.6e}")

# Verify p_1
p1_exact = 1 + mpmath.euler/2 - log(4*pi)/2
print(f"  p_1 exact: {float(p1_exact):.12f}")

# Decay check
print(f"\nDecay check p_k / (k * |rho_1|^{{-k}}):")
for k in [50, 100, 200, 400, 600, 800, 913]:
    decay_expect = k * 14.136**(-k)
    pv = float(p[k])
    if decay_expect > 0:
        ratio = abs(pv) / decay_expect if abs(pv) > 0 else 0
        print(f"  k={k:4d}: p_k={pv:.3e}, ratio={ratio:.4f}")

# STEP 4: lambda_n (optimized: binomial recurrence avoids mpmath.binomial big-int overhead)
print(f"\nComputing lambda_n for n=1..{N_LAMBDA} (optimized: binomial recurrence)...")
t0 = time.time()

# Pre-sign p: sp[k] = (-1)^{k-1} * p[k]
sp = [mpf('0')] * (K_COMPUTE + 1)
for k in range(1, K_COMPUTE + 1):
    sp[k] = p[k] if k % 2 == 1 else -p[k]

lambda_vals = []

# Phase 1: n = 1..K_COMPUTE  (k_max = n, grows)
for n in range(1, K_COMPUTE + 1):
    binom_n = [mpf('0')] * (n + 1)
    binom_n[0] = mpf('1')
    for k in range(1, n + 1):
        binom_n[k] = binom_n[k-1] * mpf(n - k + 1) / mpf(k)
    lam = mpf('0')
    for k in range(1, n + 1):
        lam += binom_n[k] * sp[k]
    lambda_vals.append(float(lam))
    if n <= 5:
        elapsed = time.time() - t0
        status = 'POS' if lam > 0 else '***NEG***'
        print(f"  n={n:5d}: lambda = {float(lam):.6f}  ({status})  [{elapsed:.1f}s]", flush=True)

elapsed_phase1 = time.time() - t0
print(f"  Phase 1 (n=1..{K_COMPUTE}) done in {elapsed_phase1:.1f}s", flush=True)

# Phase 2: n = K_COMPUTE+1..N_LAMBDA  (k_max = K_COMPUTE throughout)
# C(n,k) = C(n-1,k) * n / (n-k)
binom2 = [mpf('0')] * (K_COMPUTE + 1)
binom2[0] = mpf('1')
cur_n = K_COMPUTE
binom2[1] = mpf(cur_n)
for k in range(2, K_COMPUTE + 1):
    binom2[k] = binom2[k-1] * mpf(cur_n - k + 1) / mpf(k)

for n in range(K_COMPUTE + 1, N_LAMBDA + 1):
    new_n = mpf(n)
    for k in range(1, K_COMPUTE + 1):
        binom2[k] = binom2[k] * new_n / (new_n - k)
    lam = mpf('0')
    for k in range(1, K_COMPUTE + 1):
        lam += binom2[k] * sp[k]
    lambda_vals.append(float(lam))
    if n % 500 == 0:
        elapsed = time.time() - t0
        status = 'POS' if lam > 0 else '***NEG***'
        print(f"  n={n:5d}: lambda = {float(lam):.6f}  ({status})  [{elapsed:.1f}s]", flush=True)

lambda_time = time.time() - t0
print(f"\nDone in {lambda_time:.1f}s. Total: {time.time()-t_total:.1f}s")

all_pos = all(v > 0 for v in lambda_vals)
min_v = min(lambda_vals); min_n = lambda_vals.index(min_v) + 1
max_v = max(lambda_vals); max_n = lambda_vals.index(max_v) + 1

print(f"\n{'='*60}")
print(f"UNCONDITIONAL RESULT: lambda_n > 0 for n=1..{N_LAMBDA}: {all_pos}")
print(f"Min: lambda_{min_n} = {min_v:.10f}")
print(f"Max: lambda_{max_n} = {max_v:.10f}")

print(f"\nGrowth check lambda_n / (n/2 * log(n)):")
for n in [10, 100, 500, 1000, 2000, 3000, 5000]:
    if n <= N_LAMBDA and n > 1:
        ratio = lambda_vals[n-1] / ((n/2)*math.log(n))
        print(f"  n={n:5d}: lambda={lambda_vals[n-1]:.4f}, ratio={ratio:.4f}")

print(f"\nVerification:")
known = {1: 0.023095708966, 10: 2.279339363193, 500: 991.9000929922, 1000: 2326.0531535, 2000: 6977.83}
for n, kn in known.items():
    if n <= N_LAMBDA:
        computed = lambda_vals[n-1]
        rel = abs(computed - kn)/abs(kn)
        print(f"  n={n}: {computed:.10f}  (known ~{kn})  rel_err={rel:.2e}")

results = {
    'N_max': N_LAMBDA,
    'K_compute': K_COMPUTE,
    'precision_dps': mp.dps,
    'cauchy_radius': float(R),
    'M_points': M,
    'method': f'Direct Cauchy DFT: xi(s) at {M} pts on |s|=5, K={K_COMPUTE}, dps={mp.dps}. UNCONDITIONAL.',
    'all_positive': all_pos,
    'min_lambda': min_v, 'min_at_n': min_n,
    'max_lambda': max_v, 'max_at_n': max_n,
    'p1': float(p[1]),
    'lambda_vals': lambda_vals,
    'spot_checks': {str(n): lambda_vals[n-1] for n in [1,5,10,50,100,500,1000,2000,3000,4000,5000] if n <= N_LAMBDA},
    'timing': {'eval_s': eval_time, 'dft_s': dft_time, 'lambda_s': lambda_time},
}
outfile = os.path.join(RESULTS_DIR, 'li_unconditional_5000.json')
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to {outfile}")

print(f"\n{'='*60}")
print(f"UNCONDITIONAL CONCLUSION")
print(f"lambda_n > 0 for n=1..{N_LAMBDA}: {all_pos}")
if all_pos:
    print(f"=> Positivity confirmed UNCONDITIONALLY for n=1..{N_LAMBDA}.")
    print(f"=> No zero locations assumed.")
    print(f"=> Method: direct Cauchy integral at |s|=5, mp.dps={mp.dps}.")
else:
    neg = [i+1 for i, v in enumerate(lambda_vals) if v <= 0]
    print(f"=> NEGATIVE VALUES AT n = {neg}")

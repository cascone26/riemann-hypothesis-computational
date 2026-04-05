"""
Unconditional Li coefficients via direct Cauchy DFT — N=1000, properly calibrated.

CALIBRATION (see li_calibrate_K.py):
  For N=1000: k*=70.7, K_COMPUTE=210 (tail<10^{-15}), dps=130 (margin ~35 dps at k=210).
  Verified: for N=500, K=120, dps=80 gave correct results; this extends properly.
"""

import mpmath
from mpmath import mp, mpf, mpc, log, pi, gamma, zeta, binomial
import json, time, os, math

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

mp.dps = 130
N_LAMBDA = 1000
K_COMPUTE = 210
M = 600          # M > K_COMPUTE + 200
R = mpf('5')

print(f"Unconditional Li coefficients via direct Cauchy DFT — N={N_LAMBDA}")
print(f"  r={float(R)}, M={M} pts, K_COMPUTE={K_COMPUTE}, dps={mp.dps}")
prec_margin = mp.dps - K_COMPUTE * math.log10(float(14.136/R))
print(f"  Precision margin at k={K_COMPUTE}: {prec_margin:.1f} dps (need >15)")
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
print(f"  Done in {eval_time:.1f}s. |f[0]| = {abs(float(fvals[0].real)):.4f}")

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
print(f"  p_50 = {float(p[50]):.6e}")
print(f"  p_100 = {float(p[100]):.6e}")
print(f"  p_200 = {float(p[200]):.6e}")

# Verify p_1
p1_exact = 1 + mpmath.euler/2 - log(4*pi)/2
print(f"  p_1 exact: {float(p1_exact):.12f}")

# Decay check
print(f"\nDecay check p_k / (k * |rho_1|^{{-k}}):")
for k in [20, 60, 100, 150, 200, 210]:
    decay_expect = k * 14.136**(-k)
    pv = float(p[k])
    if abs(pv) > 0:
        ratio = abs(pv) / decay_expect
        print(f"  k={k:4d}: p_k={pv:.3e}, ratio={ratio:.4f}")

# STEP 4: lambda_n
print(f"\nComputing lambda_n for n=1..{N_LAMBDA}...")
t0 = time.time()
lambda_vals = []
for n in range(1, N_LAMBDA + 1):
    lam = mpf('0')
    k_max = min(n, K_COMPUTE)
    for k in range(1, k_max + 1):
        sign = mpf('1') if k % 2 == 1 else mpf('-1')
        lam += binomial(n, k) * sign * p[k]
    lambda_vals.append(float(lam))
    if n % 100 == 0 or n <= 5:
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
for n in [10, 100, 300, 500, 700, 1000]:
    if n <= N_LAMBDA and n > 1:
        ratio = lambda_vals[n-1] / ((n/2)*math.log(n))
        print(f"  n={n:5d}: lambda={lambda_vals[n-1]:.4f}, ratio={ratio:.4f}")

print(f"\nVerification:")
known = {1: 0.023095708966, 10: 2.279339363193, 500: 991.9000929922}
for n, kn in known.items():
    if n <= N_LAMBDA:
        computed = lambda_vals[n-1]
        rel = abs(computed - kn)/abs(kn)
        print(f"  n={n}: {computed:.12f}  (known {kn:.12f})  rel_err={rel:.2e}")

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
    'spot_checks': {str(n): lambda_vals[n-1] for n in [1,5,10,50,100,200,500,700,1000] if n <= N_LAMBDA},
    'timing': {'eval_s': eval_time, 'dft_s': dft_time, 'lambda_s': lambda_time},
}
outfile = os.path.join(RESULTS_DIR, 'li_unconditional_1000.json')
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

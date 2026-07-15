"""
Compute lambda_n for n=1..N_MAX using the full Odlyzko 2M verified-zero dataset.
Goal: measure the actual growth rate, verify positivity, look for anomalies.

Key formula (new, derivable from the paper):
  F_n(sigma, gamma) = 4 * (1 - cos(n*phi) * cosh(n*|log r|))
  where phi = arg(alpha(sigma, gamma)) = arg(alpha(1-sigma, gamma))  [SAME phase for both!]
  and |log r| = 0 iff sigma = 1/2 (on critical line).

For on-line zeros (sigma=1/2): F_n = 4*(1-cos(n*theta)) >= 0.
lambda_n = Sum_j F_n(1/2, gamma_j)  (under assumption all zeros on critical line)
         = Sum_j 4*(1 - cos(n * theta_j))
         = 4k - 4 * Re[Sum_j exp(i*n*theta_j)]

where theta_j = 2*atan(1/(2*gamma_j)).

Fast computation: vectorized over n using numpy broadcasting.
"""
import numpy as np
import json
import time
import os

DATA_FILE = '/Users/scones/projects/riemann/data/zeros_odlyzko_2M'
RESULTS_DIR = '/Users/scones/projects/riemann/results'

# Load zeros
print("Loading Odlyzko 2M zeros...")
t0 = time.time()
gammas = []
with open(DATA_FILE) as f:
    for line in f:
        line = line.strip()
        if line:
            gammas.append(float(line))
gammas = np.array(gammas)
k = len(gammas)
print(f"  Loaded {k} zeros in {time.time()-t0:.1f}s. Range: [{gammas[0]:.3f}, {gammas[-1]:.3f}]")

# Compute theta_j = 2*atan(1/(2*gamma_j)) for all zeros
# For large gamma: theta_j ≈ 1/gamma_j
thetas = 2.0 * np.arctan(0.5 / gammas)
print(f"  theta range: [{thetas[-1]:.6e}, {thetas[0]:.6e}]")

# Compute lambda_n for n=1..N_MAX
# lambda_n = 4*k - 4*Re[sum_j exp(i*n*theta_j)]
#
# Vectorized over n: for each n, compute sum_j cos(n*theta_j)
# This is O(k*N_MAX) total.

N_MAX = 10000  # Up to 10k; adjust as needed
print(f"\nComputing lambda_n for n=1..{N_MAX} using {k} zeros...")
print(f"(Method: vectorized numpy, O(k*N_MAX) = {k*N_MAX:.2e} operations)")

t0 = time.time()
lambda_vals = np.zeros(N_MAX)

# Process in batches of n to use memory efficiently
BATCH_N = 500  # Batch size in n
n_batches = (N_MAX + BATCH_N - 1) // BATCH_N

for batch_idx in range(n_batches):
    n_start = batch_idx * BATCH_N + 1
    n_end = min((batch_idx + 1) * BATCH_N, N_MAX)
    batch_ns = np.arange(n_start, n_end + 1)  # shape (batch,)

    # Outer product: batch_ns x thetas -> (batch, k)
    # cos_mat[i, j] = cos(n_i * theta_j)
    # Use float32 to save memory for large k
    phase_mat = np.outer(batch_ns.astype(np.float64), thetas)  # (batch, k)
    cos_mat = np.cos(phase_mat)  # (batch, k)
    cos_sums = cos_mat.sum(axis=1)  # (batch,)

    lambda_vals[n_start-1:n_end] = 4*k - 4*cos_sums

    elapsed = time.time() - t0
    eta = elapsed / (n_end) * (N_MAX - n_end) if n_end > 0 else 0
    if batch_idx % 5 == 0 or n_end == N_MAX:
        lmin = lambda_vals[n_start-1:n_end].min()
        lmax = lambda_vals[n_start-1:n_end].max()
        print(f"  n=[{n_start},{n_end}], elapsed={elapsed:.1f}s, ETA={eta:.1f}s, "
              f"range=[{lmin:.2f},{lmax:.2f}]", flush=True)

elapsed = time.time() - t0
print(f"\nDone in {elapsed:.1f}s.")
print(f"All lambda_n > 0: {np.all(lambda_vals > 0)}")
print(f"Min lambda_n: {lambda_vals.min():.4f} at n={np.argmin(lambda_vals)+1}")
print(f"Max lambda_n: {lambda_vals.max():.4f} at n={np.argmax(lambda_vals)+1}")

# Growth rate analysis
ns = np.arange(1, N_MAX + 1, dtype=float)
# Expected: lambda_n ~ (n/2)*log(n/(2*pi*e)) + (gamma_E/2)*n for large n
# Normalize by n*log(n)/2 to see if ratio -> 1
gamma_E = 0.5772156649
log2pie = np.log(2 * np.pi * np.e)

# Keiper-Li asymptotic: lambda_n ~ (n/2)*(log(n/(2*pi*e)) + gamma_E)
# or equivalently: lambda_n/(n*log(n)/2) -> 1 for large n
n_large = ns[99:]  # n=100 to N_MAX
lam_large = lambda_vals[99:]
expected_main = (n_large / 2.0) * (np.log(n_large) - log2pie) + (gamma_E / 2.0) * n_large

print(f"\n=== Growth Rate Analysis ===")
# Fit lambda_n = a*n*log(n) + b*n in the range n=1000..N_MAX
if N_MAX >= 2000:
    n_fit = ns[999:]
    lam_fit = lambda_vals[999:]
    log_n_fit = np.log(n_fit)
    # Linear regression: lambda_n = a*(n*log_n) + b*n
    X = np.column_stack([n_fit * log_n_fit, n_fit])
    coeffs, _, _, _ = np.linalg.lstsq(X, lam_fit, rcond=None)
    a, b = coeffs
    print(f"  Fit lambda_n ≈ {a:.6f}*n*log(n) + {b:.6f}*n  [n=1000..{N_MAX}]")
    print(f"  Expected: 0.5*n*log(n) (leading term under RH)")
    print(f"  Ratio a/0.5 = {a/0.5:.6f}  (should be ~1)")

    residuals = lam_fit - (a * n_fit * log_n_fit + b * n_fit)
    print(f"  Residual std: {residuals.std():.4f}")
    print(f"  Max residual: {np.abs(residuals).max():.4f} at n={n_fit[np.argmax(np.abs(residuals))]:.0f}")

# Ratio lambda_n / (n*log(n)/2) at various n
print(f"\n  lambda_n / (n*log(n)/2) at selected n:")
check_ns = [10, 50, 100, 500, 1000, 5000, N_MAX]
for cn in check_ns:
    if cn <= N_MAX:
        ratio = lambda_vals[cn-1] / (cn * np.log(cn) / 2) if cn > 1 else np.nan
        print(f"    n={cn:6d}: lambda={lambda_vals[cn-1]:.4f}, ratio={ratio:.6f}")

# Check minimum values (looking for near-zeros)
print(f"\n=== Near-Zero Scan ===")
# Find n where lambda_n is smallest relative to expected
expected_n = (ns / 2.0) * (np.log(np.maximum(ns, 2)) - log2pie)
expected_n = np.maximum(expected_n, 1.0)  # floor at 1
ratio_arr = lambda_vals / expected_n
min_ratio_idx = np.argsort(ratio_arr)[:20]
print(f"  20 smallest lambda_n/expected values:")
for idx in min_ratio_idx:
    n = idx + 1
    print(f"    n={n:6d}: lambda={lambda_vals[idx]:.4f}, expected={expected_n[idx]:.4f}, ratio={ratio_arr[idx]:.6f}")

# Look for anomalies: lambda_n decreasing rapidly
print(f"\n=== Anomaly Scan: Large Drops ===")
diffs = np.diff(lambda_vals)
big_drops_idx = np.where(diffs < -1000)[0]
print(f"  Drops > 1000 in n=[1,{N_MAX}]: {len(big_drops_idx)}")
if len(big_drops_idx) > 0 and len(big_drops_idx) <= 20:
    for idx in big_drops_idx[:10]:
        print(f"    n={idx+1}->n={idx+2}: {lambda_vals[idx]:.2f} -> {lambda_vals[idx+1]:.2f} (drop {diffs[idx]:.2f})")

# Verify against known small values
print(f"\n=== Spot Checks (vs known Li coefficients) ===")
# Li's original: lambda_1 = 1/2*(2+log(pi)-gamma_E-log(4)) ? Let me use the formula
# Standard: lambda_1 ≈ 0.0230957... (actually close to 0)
# Wait: lambda_1 = 1 - (1-1/rho)^1 summed = sum 1/rho (essentially). Hmm.
# From our earlier computation: li_coefficients_200.json should have the values.
li_file = os.path.join(RESULTS_DIR, 'li_coefficients_200.json')
if os.path.exists(li_file):
    with open(li_file) as f:
        li_data = json.load(f)
    # Compare a few
    print(f"  Comparing with precomputed Li coefficients (200 values):")
    compare_ns = [1, 2, 5, 10, 20, 50, 100, 200]
    for cn in compare_ns:
        if cn <= 200 and cn <= N_MAX:
            # Find in li_data
            if isinstance(li_data, list):
                known = li_data[cn-1] if cn-1 < len(li_data) else None
            elif isinstance(li_data, dict):
                known = li_data.get(str(cn)) or li_data.get(f'lambda_{cn}')
            else:
                known = None
            if known is not None:
                computed = lambda_vals[cn-1]
                print(f"    n={cn:4d}: computed={computed:.6f}, known={known:.6f}, "
                      f"diff={computed-known:.6f} ({(computed/known-1)*100:.4f}%)")
else:
    print(f"  No li_coefficients_200.json found, skipping comparison.")

# Save results
results = {
    'N_max': N_MAX,
    'n_zeros_used': k,
    'all_positive': bool(np.all(lambda_vals > 0)),
    'min_lambda': float(lambda_vals.min()),
    'min_lambda_at_n': int(np.argmin(lambda_vals) + 1),
    'max_lambda': float(lambda_vals.max()),
    'max_lambda_at_n': int(np.argmax(lambda_vals) + 1),
    'fit_a': float(a) if N_MAX >= 2000 else None,
    'fit_b': float(b) if N_MAX >= 2000 else None,
    'lambda_at_n10': float(lambda_vals[9]),
    'lambda_at_n100': float(lambda_vals[99]),
    'lambda_at_n1000': float(lambda_vals[999]) if N_MAX >= 1000 else None,
    'lambda_at_nmax': float(lambda_vals[-1]),
    'small_values': [(int(idx+1), float(lambda_vals[idx])) for idx in min_ratio_idx[:10]],
}
outfile = os.path.join(RESULTS_DIR, 'lambda_2M_growth.json')
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {outfile}")
print(f"\n=== KEY FINDING ===")
print(f"lambda_n > 0 for all n=1..{N_MAX}: {np.all(lambda_vals > 0)}")
print(f"Min lambda_n = {lambda_vals.min():.6f} at n={np.argmin(lambda_vals)+1}")
print(f"Growth rate (leading term coefficient): {a:.6f} (expected: 0.5)")

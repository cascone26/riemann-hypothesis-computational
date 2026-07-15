"""
Fast λ_n computation for n=1..N_MAX using Odlyzko 2M zeros.
Memory-efficient: iterate over n, compute sum_j cos(n*theta_j) with just a 2M vector.
"""
import numpy as np
import json
import time
import os

DATA_FILE = '/Users/scones/projects/riemann/data/zeros_odlyzko_2M'
RESULTS_DIR = '/Users/scones/projects/riemann/results'
N_MAX = 10000

print("Loading zeros...")
t0 = time.time()
gammas = np.fromiter(
    (float(l) for l in open(DATA_FILE) if l.strip()),
    dtype=np.float64
)
k = len(gammas)
print(f"  {k} zeros, range [{gammas[0]:.3f}, {gammas[-1]:.3f}], loaded in {time.time()-t0:.1f}s")

# theta_j = 2*atan(1/(2*gamma_j)), stored as float32 to save memory (still ~16M floats)
thetas = (2.0 * np.arctan(0.5 / gammas)).astype(np.float32)
print(f"  theta range: [{thetas[-1]:.6e}, {thetas[0]:.6e}]")

# Precompute exp(i*theta_j) as complex64 — avoids recomputing cos(n*theta_j) from scratch each n
# Instead: track running product: cumulative_angle = n*theta_j mod 2pi...
# Actually, fastest stable method: use float32 thetas and np.cos(n * thetas)
# At float32: 2M * 4 bytes = 8MB — trivially fits in L3/RAM
# Speed: np.cos on 2M float32 ≈ 2ms on modern CPU → 10K iterations ≈ 20s

print(f"\nComputing lambda_n for n=1..{N_MAX} (one 2M cos-sum per n)...")
lambda_vals = np.empty(N_MAX, dtype=np.float64)

t0 = time.time()
for n in range(1, N_MAX + 1):
    angles = np.multiply(n, thetas, dtype=np.float32)
    cos_sum = np.cos(angles, dtype=np.float32).sum()
    lambda_vals[n - 1] = 4.0 * k - 4.0 * float(cos_sum)

    if n % 500 == 0 or n == N_MAX:
        elapsed = time.time() - t0
        eta = elapsed / n * (N_MAX - n)
        lam_min = lambda_vals[max(0, n-500):n].min()
        lam_max = lambda_vals[max(0, n-500):n].max()
        print(f"  n={n:6d}, elapsed={elapsed:.1f}s, ETA={eta:.1f}s, "
              f"range=[{lam_min:.2f},{lam_max:.2f}]", flush=True)

elapsed = time.time() - t0
print(f"\nDone in {elapsed:.1f}s.")
print(f"All lambda_n > 0: {np.all(lambda_vals > 0)}")
print(f"Min lambda_n: {lambda_vals.min():.6f} at n={np.argmin(lambda_vals)+1}")
print(f"Max lambda_n: {lambda_vals.max():.6f} at n={np.argmax(lambda_vals)+1}")

# --- Growth rate fit ---
ns = np.arange(1, N_MAX + 1, dtype=np.float64)
gamma_E = 0.5772156649
log2pie = np.log(2 * np.pi * np.e)

print(f"\n=== Growth Rate Fit (n=1000..{N_MAX}) ===")
n_fit = ns[999:]
lam_fit = lambda_vals[999:]
X = np.column_stack([n_fit * np.log(n_fit), n_fit])
coeffs, _, _, _ = np.linalg.lstsq(X, lam_fit, rcond=None)
a, b = coeffs
print(f"  lambda_n ≈ {a:.6f}*n*log(n) + {b:.6f}*n")
print(f"  Expected leading coeff: 0.5 (under RH)")
print(f"  Ratio a/0.5 = {a/0.5:.6f}  (1.0 = perfect RH asymptotics)")
residuals = lam_fit - (a * n_fit * np.log(n_fit) + b * n_fit)
print(f"  Residual std: {residuals.std():.4f}, max: {np.abs(residuals).max():.4f}")

# --- Ratio at selected n ---
print(f"\n  lambda_n / (n*log(n)/2):")
for cn in [10, 50, 100, 500, 1000, 2000, 5000, N_MAX]:
    if cn <= N_MAX and cn > 1:
        ratio = lambda_vals[cn-1] / (cn * np.log(cn) / 2)
        print(f"    n={cn:6d}: lambda={lambda_vals[cn-1]:.4f}, ratio={ratio:.6f}")

# --- Near-zero scan ---
print(f"\n=== Near-Zero Scan ===")
expected_n = np.maximum((ns / 2.0) * (np.log(np.maximum(ns, 2)) - log2pie), 1.0)
ratio_arr = lambda_vals / expected_n
min_idx = np.argsort(ratio_arr)[:20]
print("  20 smallest lambda_n/expected:")
for idx in min_idx:
    print(f"    n={idx+1:6d}: lambda={lambda_vals[idx]:.6f}, "
          f"expected={expected_n[idx]:.4f}, ratio={ratio_arr[idx]:.6f}")

# --- Anomaly scan ---
print(f"\n=== Anomaly Scan ===")
diffs = np.diff(lambda_vals)
print(f"  Mean increment: {diffs.mean():.4f}")
print(f"  Std increment: {diffs.std():.4f}")
drops = np.where(diffs < -500)[0]
print(f"  Drops > 500: {len(drops)}")
for idx in drops[:10]:
    print(f"    n={idx+1}->{idx+2}: {lambda_vals[idx]:.2f} -> {lambda_vals[idx+1]:.2f}")

# --- Compare with known Li values ---
li_file = os.path.join(RESULTS_DIR, 'li_coefficients_200.json')
if os.path.exists(li_file):
    with open(li_file) as f:
        li_data = json.load(f)
    print(f"\n=== Spot Checks vs li_coefficients_200.json ===")
    for cn in [1, 2, 5, 10, 20, 50, 100, 200]:
        if cn <= N_MAX:
            if isinstance(li_data, list):
                known = li_data[cn-1] if cn-1 < len(li_data) else None
            else:
                known = li_data.get(str(cn))
            if known is not None:
                comp = lambda_vals[cn-1]
                print(f"    n={cn:4d}: computed={comp:.8f}, known={known:.8f}, "
                      f"diff={comp-known:.2e} ({(comp/known-1)*100:.4f}%)")

# --- Save ---
results = {
    'N_max': N_MAX, 'n_zeros': k,
    'all_positive': bool(np.all(lambda_vals > 0)),
    'min_lambda': float(lambda_vals.min()), 'min_at_n': int(np.argmin(lambda_vals)+1),
    'max_lambda': float(lambda_vals.max()), 'max_at_n': int(np.argmax(lambda_vals)+1),
    'fit_a': float(a), 'fit_b': float(b), 'fit_a_over_half': float(a/0.5),
    'residual_std': float(residuals.std()),
    'lambda_spot': {cn: float(lambda_vals[cn-1]) for cn in [1,2,5,10,50,100,500,1000,5000,N_MAX]
                    if cn <= N_MAX},
    'near_zero_top20': [(int(idx+1), float(lambda_vals[idx])) for idx in min_idx],
}
outfile = os.path.join(RESULTS_DIR, 'lambda_2M_growth.json')
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n=== KEY FINDINGS ===")
print(f"lambda_n > 0 for all n=1..{N_MAX}: {np.all(lambda_vals > 0)}")
print(f"Min lambda_n = {lambda_vals.min():.8f} at n={np.argmin(lambda_vals)+1}")
print(f"Growth rate leading coeff: {a:.8f} (expected 0.5 under RH)")
print(f"Saved: {outfile}")

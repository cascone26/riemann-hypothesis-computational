"""
Platt-Trudgian (2021) noncancellation certificate using 20M LMFDB zeros.

Formula:
  lambda_n = sum_{rho} [1 - (1-1/rho)^n]
  lambda_n^(M) = 4*K - 4*sum_j cos(n*theta_j)  where theta_j = 2*arctan(1/(2*gamma_j))

Under RH: lambda_n >= lambda_n^(M) (finite sum is a lower bound, since adding
more on-line zeros only increases lambda_n).

With Platt-Trudgian (2021) unconditional bound:
  lambda_n >= lambda_n^(M) - n^2 * c_platt
  where c_platt = C / T_M^2, T_M = height of M-th verified zero.
  C ≈ 2.147 (empirical from 2M cert check with c_platt=1.6747e-12, T_{2M}=1132490.66)

Cert: lambda_n > 0  if  lambda_n^(M) > n^2 * c_platt
"""

import numpy as np
import time
import os

DATA_DIR = '/Users/scones/projects/riemann/data'

# Files in order (sorted by zero index)
FILES = [
    (os.path.join(DATA_DIR, 'zeros_odlyzko_2M'),        'odlyzko_2M'),
    (os.path.join(DATA_DIR, 'zeros_lmfdb_2Mto20M'),     'lmfdb_2M-8M'),
    (os.path.join(DATA_DIR, 'zeros_lmfdb_8Mto14M_hp'),  'lmfdb_8M-14M'),
    (os.path.join(DATA_DIR, 'zeros_lmfdb_14Mto20M_aspire'), 'lmfdb_14M-20M'),
]

# Platt-Trudgian constant C, calibrated from 2M cert:
# c_platt(2M) = 1.6747e-12, T_{2M} = 1132490.66
# C = c_platt * T^2 = 1.6747e-12 * (1132490.66)^2 = 2.147
PLATT_C = 2.147

print("=" * 70)
print("PLATT-TRUDGIAN CERT CHECK — 20M zeros")
print("=" * 70)

# --- Load all zeros ---
print("\nLoading zeros...")
all_gammas = []
t0 = time.time()
for fpath, label in FILES:
    t1 = time.time()
    chunk = []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    val = float(line.split()[0])
                    chunk.append(val)
                except:
                    pass
    arr = np.array(chunk, dtype=np.float64)
    print(f"  {label}: {len(arr):,} zeros, range [{arr[0]:.3f}, {arr[-1]:.3f}] [{time.time()-t1:.1f}s]")
    all_gammas.append(arr)

gammas = np.concatenate(all_gammas)
K = len(gammas)
T_M = float(gammas[-1])
c_platt = PLATT_C / (T_M ** 2)
print(f"\nTotal: {K:,} zeros loaded in {time.time()-t0:.1f}s")
print(f"T_M (20M-th zero) = {T_M:,.3f}")
print(f"c_platt = {PLATT_C} / T_M^2 = {c_platt:.6e}")

# --- Compute thetas (float32 saves ~80MB RAM) ---
print("\nComputing theta angles (float32)...")
t0 = time.time()
thetas = (2.0 * np.arctan(0.5 / gammas)).astype(np.float32)
print(f"  theta range: [{float(thetas[-1]):.6e}, {float(thetas[0]):.6e}]  [{time.time()-t0:.1f}s]")
del gammas  # free RAM

# --- Cert check at sparse n values ---
print(f"\n{'n':>15}  {'lam^(20M)':>14}  {'n2*c_platt':>14}  {'ratio':>10}  cert?")
print("-" * 70)

def check_cert(n, thetas, K, c_platt):
    angles = (n * thetas.astype(np.float64)).astype(np.float32)
    cos_sum = float(np.cos(angles, dtype=np.float32).sum())
    # Correct Li formula: lambda_n = 2*K - 2*sum_j cos(n*theta_j)
    # (factor of 2: each pair rho,rho_bar contributes 2*(1-cos))
    lam = 2.0 * K - 2.0 * cos_sum
    tail = c_platt * n * n
    ratio = lam / tail if tail > 0 else float('inf')
    return lam, tail, ratio

# Sparse n grid: powers of 10, then fine-tune near the limit
import sys

sparse_ns = [
    1, 10, 100, 1_000, 10_000, 100_000, 1_000_000,
    10_000_000, 100_000_000, 1_000_000_000,
    5_000_000_000, 10_000_000_000, 20_000_000_000,
    50_000_000_000, 100_000_000_000,
]

cert_limit = None
prev_ratio = None
prev_n = None

for n in sparse_ns:
    lam, tail, ratio = check_cert(n, thetas, K, c_platt)
    cert = "YES" if ratio > 1.0 else "*** NO ***"
    print(f"{n:>15,}  {lam:>14.2f}  {tail:>14.4f}  {ratio:>10.4f}  {cert}")
    sys.stdout.flush()
    if ratio < 1.0 and cert_limit is None and prev_ratio is not None and prev_ratio > 1.0:
        cert_limit = (prev_n, n)
    prev_ratio = ratio
    prev_n = n

# --- Binary search for exact cert limit ---
if cert_limit is not None:
    lo, hi = cert_limit
    print(f"\nBinary searching cert limit in [{lo:,}, {hi:,}]...")
    while hi - lo > hi * 1e-4:
        mid = int((lo + hi) // 2)
        lam, tail, ratio = check_cert(mid, thetas, K, c_platt)
        if ratio > 1.0:
            lo = mid
        else:
            hi = mid
    print(f"\nCert limit: n ~ {lo:,} to {hi:,}")
    print(f"  At n={lo:,}: ratio = {check_cert(lo, thetas, K, c_platt)[2]:.6f}")
    print(f"  At n={hi:,}: ratio = {check_cert(hi, thetas, K, c_platt)[2]:.6f}")

print(f"\n=== RESULTS ===")
print(f"Zeros used: {K:,}")
print(f"T_M (last zero height): {T_M:,.2f}")
print(f"c_platt (20M): {c_platt:.6e}")
if cert_limit:
    print(f"CERTIFICATE: λ_n > 0 for n = 1..{lo:,} [Platt-Trudgian 2021]")
else:
    print(f"CERTIFICATE: All tested n values passed (extend grid to find limit)")

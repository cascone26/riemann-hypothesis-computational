"""
Cross-validate: Cauchy DFT (unconditional) vs Zero-location (from Odlyzko 2M zeros).

Two completely different methods:
1. Cauchy DFT: evaluates xi(s) directly on |s|=5, extracts p_k via DFT, computes lambda_n.
   - NO assumption about zero locations
   - Precision bounded by dps (150 dps for N=1000)

2. Zero-location: sum_{rho} [1 - (1 - 1/rho)^n] for first 2M zeros
   - ASSUMES zeros are on critical line (or uses actual Odlyzko positions)
   - Truncation error from missing zeros beyond 2M

If they agree to many digits: strong cross-validation of both methods.
"""

import json, math, os
import numpy as np

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

# Load Cauchy DFT results
with open(os.path.join(RESULTS_DIR, 'li_unconditional_1000.json')) as f:
    data_1000 = json.load(f)
with open(os.path.join(RESULTS_DIR, 'li_unconditional_2000.json')) as f:
    data_2000 = json.load(f)

dft_1000 = data_1000['lambda_vals']  # n=1..1000
dft_2000 = data_2000['lambda_vals']  # n=1..2000

# Load zero-location results
li_2M = np.load(os.path.join(RESULTS_DIR, 'li_100k_2M.npy'))  # n=1..100000 from 2M zeros

print("Cross-Validation: Cauchy DFT vs Zero-Location (2M Odlyzko zeros)")
print("="*70)
print()
print("Method 1 (DFT): unconditional, xi(s) on |s|=5, dps=130, no zero assumptions")
print("Method 2 (ZL):  from 2M Odlyzko zeros on critical line, float64")
print()

# Key known values for comparison
known = {
    1:   0.02309570896612103,
    2:   0.09234523552809697,
    5:   0.57554336326981406,
    10:  2.27933936319283,
    50:  43.53109608,
    100: 118.60377495,
    200: 306.65576484,
    500: 991.90009299,
    1000: 2326.0531535,
}

print(f"{'n':>6} | {'DFT (dps=130)':>18} | {'ZL (2M zeros)':>18} | {'known exact':>18} | {'DFT-known':>12} | {'ZL-known':>12}")
print("-"*90)

for n in sorted(known.keys()):
    if n <= len(dft_1000):
        dft = dft_1000[n-1]
    elif n <= len(dft_2000):
        dft = dft_2000[n-1]
    else:
        dft = float('nan')

    zl = li_2M[n-1] if n <= len(li_2M) else float('nan')
    kn = known[n]

    dft_err = abs(dft - kn) / abs(kn) if not math.isnan(dft) else float('nan')
    zl_err = abs(zl - kn) / abs(kn) if not math.isnan(zl) else float('nan')

    print(f"{n:>6d} | {dft:>18.12f} | {zl:>18.12f} | {kn:>18.12f} | {dft_err:>12.2e} | {zl_err:>12.2e}")

print()
print("Agreement between DFT and ZL methods:")
print()
print(f"{'n':>6} | {'DFT':>18} | {'ZL':>18} | {'abs diff':>14} | {'rel diff':>12}")
print("-"*75)

check_ns = [1, 2, 3, 5, 10, 20, 50, 100, 200, 300, 500, 700, 1000]
for n in check_ns:
    if n > len(dft_1000): continue
    dft = dft_1000[n-1]
    zl = li_2M[n-1]
    diff = abs(dft - zl)
    rel = diff / abs(dft) if abs(dft) > 0 else float('nan')
    print(f"{n:>6d} | {dft:>18.10f} | {zl:>18.10f} | {diff:>14.2e} | {rel:>12.2e}")

print()
print("Statistical agreement summary:")
ns_check = range(1, min(len(dft_1000), len(li_2M)) + 1)
diffs = [abs(dft_1000[n-1] - li_2M[n-1]) / abs(dft_1000[n-1]) for n in range(1, 1001)
         if abs(dft_1000[n-1]) > 1e-10]
print(f"  Relative diff DFT vs ZL over n=1..1000:")
print(f"    Median: {sorted(diffs)[len(diffs)//2]:.2e}")
print(f"    Max:    {max(diffs):.2e}")
print(f"    Min:    {min(diffs):.2e}")
print()
print("INTERPRETATION:")
print("  If rel_diff DFT vs ZL << 1e-6, the two methods are consistent.")
print("  If rel_diff ~ 1e-9..1e-11, DFT is more precise (higher dps) and ZL")
print("  is truncation-limited (missing zeros beyond 2M).")
print()
print("  Perfect agreement strengthens the unconditional claim:")
print("  if an independent zero-location method (which IS conditional on zeros")
print("  being on the line) agrees with our unconditional DFT to 10+ digits,")
print("  this cross-validates both the zero positions AND our xi evaluation.")

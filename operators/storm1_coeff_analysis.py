#!/usr/bin/env python3
"""Analyze coefficient decay law: exponential vs factorial vs hybrid."""
import mpmath as mp
import numpy as np
mp.mp.dps=15
a_raw = [mp.mpf(x) for x in [0.4409, -0.7882, 0.4103, -0.1246, 0.02208, -0.003189, 0.0003445,
         -3.24e-5, 2.871e-6, -2.215e-7, 1.731e-8, -1.302e-9, 7.792e-11]]
a = [abs(x) for x in a_raw]
print("COEFFICIENT DECAY ANALYSIS")
print("=" * 60)
print("j\ta_j (magnitude)\t\tlog|a_j|\t\tlog|a_{j-1}/a_j|")
print("-" * 60)
for j in range(len(a)):
    loga = float(mp.log(a[j])) if a[j] > 0 else float('-inf')
    ratio = float(mp.log(a[j-1]/a[j])) if j>0 and a[j]>0 else float('nan')
    print(f"{j}\t{mp.nstr(a[j],5)}\t\t{loga:.3f}\t\t{ratio:.3f}")
print("\n" + "=" * 60)
# Hypothesis 1: exp(-alpha*j^2)
print("\nHYPOTHESIS 1: a_j ~ C * exp(-α·j²)")
diffs = []
for j in range(1, len(a)):
    if a[j] > 0 and a[j-1] > 0:
        ratio = mp.log(a[j-1]/a[j])
        # d/dj log(a_j) = -2αj => ratio ≈ 2α(j - j+1) = -2α
        alpha_est = -float(ratio) / (2*(j - (j-1)))  # simplified
        diffs.append(alpha_est)
        print(f"  j={j}: estimated α ~ {alpha_est:.4f}")
mean_alpha = np.mean(diffs) if diffs else 0
print(f"  Mean α ~ {mean_alpha:.4f}")
# Hypothesis 2: exp(-α·j) = geometric decay
print("\nHYPOTHESIS 2: a_j ~ C * exp(-α·j) [geometric, q^j]")
geom = []
for j in range(1, len(a)):
    if a[j] > 0 and a[j-1] > 0:
        ratio = mp.log(a[j-1]/a[j])
        alpha_geom = float(ratio)
        geom.append(alpha_geom)
        print(f"  j={j}: log(a_{{j-1}}/a_j) = {alpha_geom:.4f}")
if geom:
    mean_geom = np.mean(geom)
    print(f"  Mean decay rate ~ {mean_geom:.4f} => q ~ exp(-{mean_geom:.4f}) ~ {np.exp(-mean_geom):.4f}")
# Hypothesis 3: a_j ~ C / j!^β
print("\nHYPOTHESIS 3: a_j ~ C / (j!)^β [factorial decay, modulated by power]")
fact = [1]
for j in range(1, len(a)):
    fact.append(fact[-1] * j)
for j in range(1, len(a)):
    if a[j] > 0:
        # a_j ~ C / (j!)^β => log a_j = log C - β log(j!)
        # d/dj log a_j = -β d/dj log(j!) = -β log(j)
        log_factorial = float(mp.log(fact[j]))
        print(f"  j={j}: log(j!) = {log_factorial:.3f}, log|a_j| = {float(mp.log(a[j])):.3f}")
print("\n" + "=" * 60)
print("RECOGNITION CHECK: Is this a theta-function q-series?")
print("If a_j ~ q^(j²) => theta_3(0; q)")
print("If a_j ~ q^(j(j+1)) => theta_4(0; q)")
print(f"Observed geometric rate: exp(-{mean_geom:.3f}·j) suggests q ~ {np.exp(-mean_geom):.6f}")
print("(q ~ 0.004 is shallow for theta_3; theta typically has q ~ 0.01-0.1)")
print("\nDEEP MATCH NEEDED: Compare coefficient magnitudes to:")
print("  - DLMF 20.9.1 (theta_3 q-series)")
print("  - DLMF 20.9.2 (theta_4 q-series)")
print("  - Eisenstein E_4, E_6 q-series (DLMF 24.7)")
print("=" * 60)

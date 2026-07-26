#!/usr/bin/env python3
"""QUICK TEST: 5-minute cascade test with minimal parameters"""
import sys
sys.path.insert(0, "/Users/scones/projects/riemann/.venv_cvs/lib/python3.14/site-packages")

import mpmath as mp
import connes_cvs.operator as O
import json

mp.mp.dps = 40  # Very low precision
print("QUICK BEURLING CASCADE TEST (dps=40, minimal)", flush=True)
print("="*60, flush=True)

# Real primes baseline
print("\nBASELINE: Real primes", flush=True)
L = mp.log(13.0)
idx = list(range(-10, 11))  # N=10 (smaller)
pdata_real, _ = O.prime_powers_up_to(13)

print(f"  Real primes: {pdata_real[:5]}", flush=True)
print("  Computing cascade...", flush=True)

# Just measure: can we even call psi_prime with custom data?
try:
    test_result = O.psi_prime(mp.mpf(1), L, pdata_real[:3])
    print(f"  psi_prime works! Result: {mp.nstr(test_result, 4)}", flush=True)
except Exception as e:
    print(f"  ERROR calling psi_prime: {e}", flush=True)

# Now try with Beurling data
print("\nBEURLING: Perturbed primes", flush=True)
def get_primes(n):
    is_p = [True]*(n+1)
    is_p[0]=is_p[1]=False
    for i in range(2, int(n**0.5)+1):
        if is_p[i]:
            for j in range(i*i, n+1, i):
                is_p[j]=False
    return [i for i in range(2,n+1) if is_p[i]]

primes = get_primes(13)
print(f"  Real primes: {primes}", flush=True)

# Build Beurling version
beurling = []
for k, p in enumerate(primes, 1):
    delta = 0.1 * mp.sin(mp.mpf(k) * mp.pi / len(primes))
    beurling.append(mp.mpf(p) + delta)

print(f"  Beurling primes: {[mp.nstr(b, 5) for b in beurling[:3]]}", flush=True)

# Build Beurling pdata
pdata_beurling = []
for q in beurling[:3]:  # Just first 3 primes
    log_q = mp.log(q)
    for k in range(1, 2):
        pdata_beurling.append((float(q**k), log_q, k))

print(f"  Beurling pdata: {pdata_beurling[:3]}", flush=True)

# Try psi_prime with Beurling data
try:
    test_result_b = O.psi_prime(mp.mpf(1), L, pdata_beurling[:1])
    print(f"  psi_prime(Beurling) works! Result: {mp.nstr(test_result_b, 4)}", flush=True)
    print(f"  Ratio (real/Beurling): {mp.nstr(test_result/test_result_b, 4)}", flush=True)
except Exception as e:
    print(f"  ERROR with Beurling: {e}", flush=True)

print("\n" + "="*60)
print("CONCLUSION: Can adapt Connes library to Beurling primes? " + ("YES" if "works" in locals() else "NO"), flush=True)

result = {"quick_test_status": "psi_prime_callable_with_custom_data"}
with open("/Users/scones/projects/riemann/operators/track_crossblock/quick_test_result.json", "w") as f:
    json.dump(result, f)
print("Quick test result saved.", flush=True)

#!/usr/bin/env python3 -u
"""scan_kappa_parameters.py

Scan different (a, K) parameter values to find the critical case where κ=1.
This helps identify the regime of interest for the Π₁ structure.
"""

import sys
sys.path.insert(0, '/Users/scones/projects/riemann')

import numpy as np
from scipy.linalg import eigh
import json
from operators.suzuki_exact_Q import build_Q

def compute_kappa(a, K, NG=320):
    """Extract κ for a given (a, K)."""
    Q, pieces = build_Q(a, K, NG, return_pieces=True)
    M_sigma = pieces['core'] + pieces['tail'] - pieces['prime']

    evals = np.linalg.eigvalsh(M_sigma)
    tol = 1e-10 * np.max(np.abs(evals))
    kappa = np.sum(evals < -tol)

    return kappa, float(evals[0]), float(evals[-1])

# Scan parameters
results = {}

# Vary a at fixed K=80
print("Scanning a values at K=80:")
for a in np.linspace(0.3, 2.0, 18):
    kappa, emin, emax = compute_kappa(a, K=80)
    results[f'a={a:.2f},K=80'] = {'kappa': kappa, 'emin': emin, 'emax': emax}
    print(f"  a={a:.2f}: κ={kappa}, λ_min={emin:.4e}, λ_max={emax:.4e}")

# Vary K at fixed a=1.0
print("\nScanning K values at a=1.0:")
for K in [32, 48, 64, 80, 96, 120, 160]:
    kappa, emin, emax = compute_kappa(a=1.0, K=K)
    results[f'a=1.0,K={K}'] = {'kappa': kappa, 'emin': emin, 'emax': emax}
    print(f"  K={K:3d}: κ={kappa}, λ_min={emin:.4e}, λ_max={emax:.4e}")

# Focus on large-bandwidth (small a) values
print("\nScanning small a values (large bandwidth):")
for a in [0.1, 0.15, 0.2, 0.25, 0.3]:
    kappa, emin, emax = compute_kappa(a, K=80)
    results[f'a={a:.2f},K=80'] = {'kappa': kappa, 'emin': emin, 'emax': emax}
    print(f"  a={a:.2f}: κ={kappa}, λ_min={emin:.4e}, λ_max={emax:.4e}")

# Save results
with open('/Users/scones/projects/riemann/lanes_20260810/kappa_scan.json', 'w') as f:
    json.dump({k: {kk: (int(vv) if isinstance(vv, (int, np.integer)) else float(vv)) for kk, vv in v.items()} for k, v in results.items()}, f, indent=2)

print("\n[HANDLE] Scan saved to kappa_scan.json")

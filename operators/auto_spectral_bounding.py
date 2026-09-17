#!/usr/bin/env python3
import sys; sys.path.insert(0, "/home/cobo/projects/riemann/operators")
from critical_eigvec_id import build_tau_np
import numpy as np

N = 64
c_values = [8, 12, 14]
results = {}

def spectral_bound(tau):
    eigvals, _ = np.linalg.eigh(tau)
    return np.min(eigvals), np.max(eigvals)

for c in c_values:
    tau, ns = build_tau_np(N, c)
    min_eig, max_eig = spectral_bound(tau)
    results[c] = (min_eig, max_eig)

print(results)
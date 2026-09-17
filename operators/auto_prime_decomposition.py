#!/usr/bin/env python3
import sys; sys.path.insert(0, "/home/cobo/projects/riemann/operators")
from critical_eigvec_id import build_tau_np
import numpy as np
from scipy.linalg import sqrtm

def construct_prime_decomposition(tau):
    M = sqrtm(tau)
    return M @ M.T  # Placeholder for actual prime-based construction

tau_values = {}
for c in [8, 10, 12, 14]:
    tau, ns = build_tau_np(64, c)
    decomposition = construct_prime_decomposition(tau)
    min_eigenvalue = np.min(np.linalg.eigvalsh(decomposition))
    tau_values[c] = {'min_eigenvalue': min_eigenvalue}

print(tau_values)
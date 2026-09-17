#!/usr/bin/env python3
import sys; sys.path.insert(0, "/home/cobo/projects/riemann/operators")
from critical_eigvec_id import build_tau_np
import numpy as np

N = 64
c_values = [8, 12, 14]
results = {}

for c in c_values:
    tau, ns = build_tau_np(N, c)
    H = np.linalg.qr(tau)[0]  # QR decomposition to find a possible extension matrix
    M = H.T @ tau @ H
    min_eig_M = np.min(np.real(np.linalg.eigvals(M)))
    results[c] = min_eig_M

print(results)
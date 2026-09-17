#!/usr/bin/env python3
import sys; sys.path.insert(0, "/home/cobo/projects/riemann/operators")
from critical_eigvec_id import build_tau_np
import numpy as np

N = 64
c_values = [8, 12, 14]
results = {}

def compute_cohomology_class(tau):
    eigvals, eigvecs = np.linalg.eigh(tau)
    cohomology_class = eigvecs[:, eigvals > 0] @ np.diag(np.sqrt(eigvals[eigvals > 0]))
    return cohomology_class

for c in c_values:
    tau, ns = build_tau_np(N, c)
    cohomology_class = compute_cohomology_class(tau)
    M = cohomology_class.T @ tau @ cohomology_class
    min_eig_M = np.min(np.real(np.linalg.eigvals(M)))
    results[c] = min_eig_M

print(results)
#!/usr/bin/env python3
import sys; sys.path.insert(0, "/home/cobo/projects/riemann/operators")
from critical_eigvec_id import build_tau_np
import numpy as np

def compute_cohomological_signature(tau):
    # Placeholder for actual cohomological computation
    return np.linalg.det(tau)

test_values = [8, 12, 14]
results = []
for c in test_values:
    tau, ns = build_tau_np(64, c)
    signature = compute_cohomological_signature(tau)
    results.append((c, signature))

print(results)
#!/usr/bin/env python3
import sys; sys.path.insert(0, "/home/cobo/projects/riemann/operators")
from critical_eigvec_id import build_tau_np
import numpy as np
from scipy.linalg import eigvalsh

def compute_intersection_signature(tau):
    eigenvalues = eigvalsh(tau)
    positive_count = np.sum(eigenvalues > 0)
    negative_count = np.sum(eigenvalues < 0)
    return {'positive': positive_count, 'negative': negative_count}

tau_values = {}
for c in [8, 10, 12, 14]:
    tau, ns = build_tau_np(64, c)
    intersection_sig = compute_intersection_signature(tau)
    tau_values[c] = {'intersection_signature': intersection_sig}

print(tau_values)
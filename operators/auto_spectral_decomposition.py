#!/usr/bin/env python3
import sys; sys.path.insert(0, "/home/cobo/projects/riemann/operators")
from critical_eigvec_id import build_tau_np
import numpy as np
tau, ns = build_tau_np(64, 14)
eigenvalues, eigenvectors = np.linalg.eigh(tau)
decomposition_error = np.linalg.norm(np.dot(eigenvectors * eigenvalues, eigenvectors.T) - tau)
print(f"Decomposition error: {decomposition_error}")
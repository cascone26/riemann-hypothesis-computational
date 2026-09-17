#!/usr/bin/env python3
import sys; sys.path.insert(0, "/home/cobo/projects/riemann/operators")
from critical_eigvec_id import build_tau_np
import numpy as np
from scipy.spatial import ConvexHull

N = 64
c_values = [8, 12, 14]
results = {}

def construct_lattice(tau):
    points = np.random.randn(2*N+1, N)
    hull = ConvexHull(points)
    lattice_basis = points[hull.vertices]
    return lattice_basis

for c in c_values:
    tau, ns = build_tau_np(N, c)
    lattice_basis = construct_lattice(tau)
    M = lattice_basis.T @ np.diag(tau) @ lattice_basis
    min_eig_M = np.min(np.real(np.linalg.eigvals(M)))
    results[c] = min_eig_M

print(results)
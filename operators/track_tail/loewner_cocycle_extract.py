#!/usr/bin/env python3
"""Verify the Fourier-Loewner form of the Q tail (see tail_QFourier_loewner_note.md).
R[k,l] = Q[k,l](k^2-l^2)/(kl) must satisfy the cocycle R[k,l]+R[l,m]+R[m,k]=0 for all triples
iff Q[k,l] = c*kl*(phi(l)-phi(k))/(k^2-l^2). Then extracts phi and the tail sequence m*psi(m)."""
import numpy as np
Q = np.load('/Users/scones/projects/riemann/results/Q51_a0549_nl4000_float.npy')
modes = list(range(1, 52, 2)); n = len(modes)
R = np.zeros((n, n))
for i, k in enumerate(modes):
    for j, l in enumerate(modes):
        if i != j: R[i, j] = Q[i, j]*(k*k - l*l)/(k*l)
worst = max(abs(R[i,j]+R[j,m]+R[m,i]) for i in range(n) for j in range(i+1,n) for m in range(j+1,n))
print(f"cocycle max violation over all triples: {worst:.3e}")
phi = np.array([R[0, j] if j else 0.0 for j in range(n)])
tail_mean = phi[-6:].mean()
print(f"phi_inf (m=41..51 mean): {tail_mean:+.5f}")
for j, m in enumerate(modes):
    if m >= 15: print(f"  m={m:2d}  phi={phi[j]:+.5f}  m*psi={m*(phi[j]-tail_mean):+.3f}")

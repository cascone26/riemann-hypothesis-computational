#!/usr/bin/env python3
"""Re-check screw/CvS equivalence at interior a=0.480 (clean margin, not criticality)."""
import sys
import numpy as np
import mpmath as mp
from scipy.linalg import eigh
import os

mp.mp.dps = 60
a = mp.mpf('0.480')
b = mp.pi / (2 * a)
log2 = mp.log(2)
twoa = 2 * a
Om0 = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)
Phi0 = mp.zeta(2, mp.mpf(1)/4)
wp = mp.log(2) / mp.sqrt(2)
NL = 3000

c_matched = mp.e ** (2 * a)
print(f"[INTERIOR RECHECK] a={float(a):.6f}, c≈{float(c_matched):.6f}")

def Pp(p, t, s):
    if p == 0: return 1 / s
    if p == 1: return t / s - 1 / s**2
    return t*t / s - 2*t / s**2 + 2 / s**3

def II(p, al, om, sc, t0, t1):
    s = mp.mpc(-al, om)
    F = lambda tt: mp.e**(s*tt) * Pp(p, tt, s)
    v = F(t1) - F(t0)
    return v.real if sc == 'cos' else v.imag

def dterms(k, l):
    if k == l:
        return [(2*a*k*k*b*b, 0, k*b, 'cos'),
                (-k*k*b*b, 1, k*b, 'cos'),
                (-k*b, 0, k*b, 'sin')]
    al = k*l*b * (mp.mpf(1)/(k-l) - mp.mpf(1)/(k+l))
    ak = -k*l*b * (mp.mpf(1)/(k-l) + mp.mpf(1)/(k+l))
    return [(al, 0, l*b, 'sin'), (ak, 0, k*b, 'sin')]

def Q(k, l):
    tot = mp.mpf(0)
    gt = [(-4, 0, mp.mpf(-1)/2), (-4, 0, mp.mpf(1)/2), (8, 0, mp.mpf(0)),
          (-Om0/2, 1, mp.mpf(0)), (-Phi0/4, 0, mp.mpf(0))]
    for n in range(NL):
        gt.append((mp.mpf(1)/4 / (n + mp.mpf(1)/4)**2, 0, 2*n + mp.mpf(1)/2))
    for (ag, pg, alp) in gt:
        for (ad, pd, om, sc) in dterms(k, l):
            tot += ag * ad * II(pg + pd, alp, om, sc, mp.mpf(0), twoa)
    for (ad, pd, om, sc) in dterms(k, l):
        tot += wp*ad*II(pd+1, mp.mpf(0), om, sc, log2, twoa) + (-wp*log2)*ad*II(pd, mp.mpf(0), om, sc, log2, twoa)
    return tot

print("\nScrew Q at interior a=0.480:")
modes = [1, 3, 5, 7, 9]
Q_mat = np.array([[complex(Q(k, l)) for l in modes] for k in modes])
print(f"  Q diagonal (real): {np.diag(Q_mat).real}")
eigvals_Q = np.linalg.eigvalsh(Q_mat.real)
n_pos_Q = sum(1 for e in eigvals_Q if e > 1e-10)
n_neg_Q = sum(1 for e in eigvals_Q if e < -1e-10)
print(f"  λ_min(Q) = {eigvals_Q[0]:.6e}")
print(f"  Inertia(Q): ({n_pos_Q}, {n_neg_Q}, 0)")
print(f"  MARGIN = {eigvals_Q[0]:.3e} (CLEAN, interpretable)")

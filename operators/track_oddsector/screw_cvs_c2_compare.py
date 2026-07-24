#!/usr/bin/env python3
"""Compare screw Q vs Groskin τ at c=2 (interior point: a≈0.3466)."""
import sys, os, numpy as np, mpmath as mp
sys.path.insert(0, 'rigor')
from screw_cvs_equiv import build_groskin_tau_mpmath
mp.mp.dps = 60

# Interior point: c=2 ⟹ a=log(2)/2
c_int = 2
a_int = mp.log(mp.mpf(c_int)) / 2
print(f"[C=2 INTERIOR COMPARE] c={c_int}, a≈{float(a_int):.6f}")

# Screw Q(k,l) at a=log(2)/2
b = mp.pi / (2 * a_int)
log2 = mp.log(2)
twoa = 2 * a_int
Om0 = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)
Phi0 = mp.zeta(2, mp.mpf(1)/4)
wp = log2 / mp.sqrt(2)
NL = 2000

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
        return [(2*a_int*k*k*b*b, 0, k*b, 'cos'), (-k*k*b*b, 1, k*b, 'cos'), (-k*b, 0, k*b, 'sin')]
    al = k*l*b * (mp.mpf(1)/(k-l) - mp.mpf(1)/(k+l))
    ak = -k*l*b * (mp.mpf(1)/(k-l) + mp.mpf(1)/(k+l))
    return [(al, 0, l*b, 'sin'), (ak, 0, k*b, 'sin')]

def Q(k, l):
    tot = mp.mpf(0)
    gt = [(-4, 0, mp.mpf(-1)/2), (-4, 0, mp.mpf(1)/2), (8, 0, mp.mpf(0)), (-Om0/2, 1, mp.mpf(0)), (-Phi0/4, 0, mp.mpf(0))]
    for n in range(NL):
        gt.append((mp.mpf(1)/4 / (n + mp.mpf(1)/4)**2, 0, 2*n + mp.mpf(1)/2))
    for (ag, pg, alp) in gt:
        for (ad, pd, om, sc) in dterms(k, l):
            tot += ag * ad * II(pg + pd, alp, om, sc, mp.mpf(0), twoa)
    for (ad, pd, om, sc) in dterms(k, l):
        tot += wp*ad*II(pd+1, mp.mpf(0), om, sc, log2, twoa) + (-wp*log2)*ad*II(pd, mp.mpf(0), om, sc, log2, twoa)
    return tot

print("\nScrew Q at c=2 (a≈0.3466):")
modes = [1, 3, 5]
Q_mat = np.array([[complex(Q(k, l)) for l in modes] for k in modes])
print(f"  Modes: {modes}")
print(f"  Q diagonal (real): {np.diag(Q_mat).real}")
eigvals_Q = np.linalg.eigvalsh(Q_mat.real)
n_pos_Q = sum(1 for e in eigvals_Q if e > 1e-10)
n_neg_Q = sum(1 for e in eigvals_Q if e < -1e-10)
print(f"  λ_min(Q) = {eigvals_Q[0]:.6e}")
print(f"  Inertia(Q): ({n_pos_Q}, {n_neg_Q}, 0)")

print("\nGroskin τ at c=2, N=2:")
tau = build_groskin_tau_mpmath(mp.mpf(c_int), 2)
print(f"  τ shape: {tau.shape}")
eigvals_tau = np.linalg.eigvalsh(tau[:3, :3].astype(float))  # 3×3 to match modes
n_pos_tau = sum(1 for e in eigvals_tau if e > 1e-10)
n_neg_tau = sum(1 for e in eigvals_tau if e < -1e-10)
print(f"  τ diagonal (real) [3×3]: {np.diag(tau[:3, :3]).astype(float)}")
print(f"  λ_min(τ) = {eigvals_tau[0]:.6e}")
print(f"  Inertia(τ): ({n_pos_tau}, {n_neg_tau}, 0)")

print("\n[RESULT]")
if n_pos_Q == n_pos_tau and n_neg_Q == n_neg_tau:
    print("  ✓ INERTIAS MATCH")
    if abs(eigvals_Q[0] - eigvals_tau[0]) / max(abs(eigvals_Q[0]), abs(eigvals_tau[0])) < 0.1:
        print("  ✓ λ_min AGREE (within 10%)")
        print("  CONCLUSION: Same form, different basis (EQUIVALENT)")
    else:
        print(f"  ✗ λ_min differ significantly: Q={eigvals_Q[0]:.3e} vs τ={eigvals_tau[0]:.3e}")
else:
    print(f"  ✗ INERTIAS DIFFER: Q{(n_pos_Q, n_neg_Q)} vs τ{(n_pos_tau, n_neg_tau)}")
    print("  CONCLUSION: Different forms OR basis effects (Gap 3 UNRESOLVED, needs deeper analysis)")

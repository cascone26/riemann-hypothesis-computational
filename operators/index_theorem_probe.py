#!/usr/bin/env python3
"""§946 THE INDEX FRAMING (roundtable §945 survivor: encode RH as a ROBUST INTEGER, not a fragile margin).

The §944 fragility finding says the Weil MARGIN (min eigenvalue) is knife-edge — generically negative,
razor-positive only at exact prime positions, breaks at delta=0.001. A margin is not a good invariant.
BUT the INERTIA of a symmetric matrix (its count of negative eigenvalues, the Morse index) is an INTEGER
that is ROBUST: it does NOT change under small perturbations unless an eigenvalue actually crosses zero.

CLAIM to test (function-field case, where RH is PROVEN, so we KNOW the ground truth):
    #neg( T_N )  =  (number of Frobenius eigenvalues OFF the circle)   [once N >= 2g, robust in N/radius/angle]
i.e. the inertia of the Weil-Toeplitz form is a topological index that COUNTS RH violations.
RH_ff  <=>  index(T) = 0.  If this integer is robust (invariant to how far off, to N, to angles), it is the
"protected integer" the §945 filter selected — the razor's fragility lives only in the LAST eigenvalue to
cross, but the COUNT is stable.

This is the clean testbed BEFORE claiming anything for number fields: verify the index framing where the
answer is known. HANDLE-tier (my own code, reproducible). No fabrication — print raw inertia, let it speak.
"""
import numpy as np
np.set_printoptions(precision=4, suppress=True)

def toeplitz(radii, angles, q, N):
    """Hermitian Weil-Toeplitz T[m,n]=s_{m-n}, s_k = normalized power sum of Frobenius eigenvalues."""
    def s(k):
        return sum((r/np.sqrt(q))**abs(k)*np.exp(1j*k*t) for r, t in zip(radii, angles))
    T = np.array([[s(m-n) for n in range(N)] for m in range(N)])
    return 0.5*(T + T.conj().T)

def inertia(T, tol=1e-9):
    """(#neg, #zero, #pos) eigenvalues, tol relative to spectral scale."""
    ev = np.linalg.eigvalsh(T)
    scale = max(1.0, np.abs(ev).max())
    neg = int(np.sum(ev < -tol*scale))
    pos = int(np.sum(ev >  tol*scale))
    zero = len(ev) - neg - pos
    return neg, zero, pos, ev

q = 25; g = 4                         # genus-4 curve / F_25: 2g = 8 Frobenius eigenvalues on |a|=5
base_ang = [0.4, 1.1, 1.9, 2.7]
ang8 = base_ang + [-a for a in base_ang]   # conjugate-symmetric, 8 eigenvalues
R = np.sqrt(q)                        # on-circle radius

print("="*94)
print(f"§946 INDEX FRAMING — function-field Weil form (curve/F_{q}, genus {g}, 2g={2*g}). RH_ff proven ground truth.")
print("Test: is #neg(T_N) a ROBUST INTEGER = number of eigenvalues moved off the circle?")
print("="*94)

# --- Test A: move k eigenvalues off-circle, see #neg vs k, across N -------------------------------
print("\n[A] #neg(T_N) as we move k of the 8 eigenvalues OFF the circle (radius 5.6). Columns = N.")
Ns = [8, 9, 10, 12, 14, 16]
print(f"{'k off-circle':>13} | " + " ".join(f"N={n:>2}" for n in Ns))
for k in range(0, 5):
    row = []
    for N in Ns:
        radii = [R]*8
        for i in range(k):
            radii[i] = 5.6        # push eigenvalue i off the circle (and its conjugate stays? no—move independently)
        neg, z, p, _ = inertia(toeplitz(radii, ang8, q, N))
        row.append(neg)
    print(f"{k:>13} | " + " ".join(f"{v:>4}" for v in row))

# --- Test B: robustness of the index to HOW FAR off-circle (should NOT depend on distance) --------
print("\n[B] Robustness: move exactly k=2 eigenvalues off-circle at varying radius. #neg should be CONSTANT.")
print(f"{'radius':>8} | {'min eig (margin, FRAGILE)':>26} | {'#neg (index, ROBUST?)':>22}")
for rad in [5.001, 5.01, 5.1, 5.6, 7.0, 12.0]:
    radii = [R]*8; radii[0] = rad; radii[1] = rad
    neg, z, p, ev = inertia(toeplitz(radii, ang8, q, 12))
    print(f"{rad:>8.3f} | {ev.min():>26.3e} | {neg:>22}")

# --- Test C: robustness to ANGLES and to random on-circle config (control: k=0 => index 0) --------
print("\n[C] Control — ALL on-circle (RH holds), random angles. #neg must be 0 (the margin may be ~0, fragile).")
rng = np.random.default_rng(7)
for trial in range(5):
    a4 = np.sort(rng.uniform(0.2, 3.0, 4))
    a8 = list(a4) + [-a for a in a4]
    neg, z, p, ev = inertia(toeplitz([R]*8, a8, q, 12))
    print(f"  trial {trial}: angles={np.round(a4,2)}  min eig={ev.min():>+.3e}  #neg={neg}  #zero={z}")

print("="*94)
print("READ: if column-A #neg locks to a k-only integer (robust across N once N>=2g) and [B] holds it constant")
print("across radius while the margin swings, the inertia IS the protected integer. [C] must give #neg=0 (RH).")
print("This is the function-field ground-truth check of the §945 'encode RH as a robust index' direction.")
print("="*94)

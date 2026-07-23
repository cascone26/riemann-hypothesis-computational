#!/usr/bin/env python3
"""§178e integrity check: is the k_λ construction sound at every c, or does the argmax-parity mode-selection
misfire (→ wrong h_0/h_4 → bogus dist, explaining the non-monotone dip at c=6)? Pure numpy, seconds.
For each c: rebuild the prolate operator, report the modes the heuristic picks (i0,i4), their TRUE parity
(even ⟺ only even-degree Legendre coeffs nonzero), and whether the vanishing-integral condition really holds."""
import numpy as np
from numpy.polynomial.legendre import legval

def _y2(Nleg):
    Y = np.zeros((Nleg + 2, Nleg + 2))
    for n in range(Nleg + 1):
        if n + 1 <= Nleg + 1: Y[n + 1, n] = (n + 1) / (2 * n + 1) * np.sqrt((2 * n + 1) / (2 * n + 3))
        if n - 1 >= 0:        Y[n - 1, n] = n / (2 * n + 1) * np.sqrt((2 * n + 1) / (2 * n - 1))
    return (Y[:Nleg, :Nleg]) @ (Y[:Nleg, :Nleg])

def true_parity(col):
    """Even iff odd-degree Legendre coeffs are ~0; odd iff even-degree ~0. Returns +1 even, -1 odd, 0 mixed."""
    ev = np.sum(np.abs(col[0::2])); od = np.sum(np.abs(col[1::2]))
    if od < 1e-6 * ev: return +1
    if ev < 1e-6 * od: return -1
    return 0

for c in [4, 5, 6, 7, 8, 9, 10]:
    lam = np.sqrt(c)
    Nleg = max(200, int(6 * 2 * np.pi * lam ** 2)); cband = 2 * np.pi * lam ** 2; n = np.arange(Nleg)
    A = np.diag((n * (n + 1)).astype(float)) + cband ** 2 * _y2(Nleg); A = (A + A.T) / 2
    w, V = np.linalg.eigh(A)
    # heuristic parity (as used in k_lambda_on) vs TRUE parity
    par_h = np.array([(-1) ** int(np.argmax(np.abs(V[:, k]))) for k in range(V.shape[1])])
    ei_h = np.where(par_h > 0)[0]
    i0h, i4h = ei_h[0], ei_h[2]
    par_t = np.array([true_parity(V[:, k]) for k in range(min(12, V.shape[1]))])
    ei_t = [k for k in range(len(par_t)) if par_t[k] > 0]
    i0t, i4t = ei_t[0], ei_t[2]
    # vanishing-integral combo as built in k_lambda_on (using heuristic modes)
    b0 = (-V[0, i0h] / V[0, i4h]) if abs(V[0, i4h]) > 1e-30 else 0.0
    # integral of the combo over [-1,1] = 2*(0th Legendre coeff) → should be ~0 by construction
    def coeff(colidx): col = V[:, colidx]; return np.array([col[m] * np.sqrt((2 * m + 1) / 2) for m in range(Nleg)])
    c0 = coeff(i0h); c4 = coeff(i4h)
    integ = 2.0 * (1.0 * c0[0] + b0 * c4[0])         # ∫ combo dy on [-1,1]
    match = "OK" if (i0h, i4h) == (i0t, i4t) else f"MISMATCH heuristic({i0h},{i4h}) vs true({i0t},{i4t})"
    print(f"c={c:>2}  heuristic even-modes i0,i4=({i0h},{i4h})  true=({i0t},{i4t})  {match}")
    print(f"      first even-mode indices (true): {ei_t[:4]}   ∫(vanishing-int combo)={integ:.2e}  b0={b0:.4f}")

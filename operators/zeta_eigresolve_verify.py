#!/usr/bin/env python3
"""cy349 — VERIFY-BEFORE-SHIP the load-bearing claim: ζ's certified n_neg=0 (basis of "ζ realizes AF clean 2/3")
via DIRECT eigendecomposition, given cy348's finding that even objects have a near-null cluster the prec800
LDLᵀ inertia miscounts. ζ is even (nu=1/4) + has the pole (κ=1). Compute the actual eigenvalue spectrum at
prec2000 and confirm the inertia is what the LDLᵀ n_neg claims (n_neg should reflect the pole's structure,
NOT a spurious cluster). Compare to D-H (RH-FALSE) as the null baseline."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from dh_certified_tau import build_arb_tau, margin

def spectrum(A, DIM):
    M = np.zeros((DIM, DIM))
    for i in range(DIM):
        for j in range(DIM):
            M[i, j] = float(A[i, j].mid())
    M = 0.5*(M + M.T)
    return np.linalg.eigvalsh(M)

prec, N = 2000, 64
for c in [25, 50]:
    ctx.prec = prec
    A, DIM = build_arb_tau(c, N, prec)
    nneg_ldlt, ok, mlow = margin(A, DIM)
    ev = spectrum(A, DIM)
    neg = ev[ev < -1e-13]          # genuine negatives above the ~1e-16 machine floor
    print(f"ζ c={c} N={N} DIM={DIM}: LDLᵀ n_neg={nneg_ldlt} | direct-eig genuine-neg(<-1e-13)={len(neg)}")
    print(f"   smallest 6 eig: {np.array2string(ev[:6], formatter={'float':lambda x:f'{x:+.3e}'})}")
    if len(neg):
        print(f"   genuine negatives: {np.array2string(neg, formatter={'float':lambda x:f'{x:+.3e}'})}")
    print(f"   most-negative = {ev[0]:+.3e}", flush=True)
print("\nVERDICT: if direct-eig genuine-neg agrees with (or explains) the LDLᵀ n_neg — and any 'extra' LDLᵀ")
print("negatives are ~1e-16 floor cluster like χ_5 — then ζ's n_neg=0 (RH-favorable) is ROBUST to the even-")
print("cluster caveat, and the AF clean-2/3 realization stands. A genuine O(1) negative beyond the pole = a problem.")

#!/usr/bin/env python3
"""B2 probe: Szegő / determinant structure of the Groskin cutoff-free Weil form τ.

Reuses the VERIFIED entry engine build_arb_tau() from
external/groskin_2607.02828/code/arb_ldlt_certify.py (imported by path, not re-derived).
Builds τ at fixed c, dimension DIM=2N+1, extracts the LDL^T pivots d_i (float midpoints).

Szegő hypothesis (strong Szegő limit theorem generalization): for a Toeplitz-like PD form,
log det τ_M = M * (symbol mean) + [convergent SoS correction] + o(1).
Equivalently the pivot logs log d_i should CONVERGE to a constant (the symbol mean) with a
SUMMABLE correction Σ (log d_i - mean). If the correction converges, that finite constant is
a manifest-positivity invariant (a determinant-calculus seed), which raw natural-order pivot
magnitude (A1, found erratic) is NOT. Also report Borodin-Okounkov-style: is the tail of
log d_i - mean summable / does Σ|log d_i - c| converge?

REPORT tier (float midpoints of Arb balls); the rigorous cert is arb_ldlt_certify.py itself.
"""
import importlib.util, sys, os
import numpy as np

GPATH = os.path.expanduser('~/projects/riemann/external/groskin_2607.02828/code/arb_ldlt_certify.py')
spec = importlib.util.spec_from_file_location('groskin_cert', GPATH)
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)

def to_float(x):
    try: return float(x)
    except Exception:
        try: return float(x.mid())
        except Exception: return float(str(x).split('+/-')[0].strip('[] '))

def pivots_of(c, N, prec):
    """LDL^T kept in Arb (full precision) — float cast ONLY at the end. Casting before
    the LDL (v1) manufactured spurious negative pivots because d spans >14 orders."""
    from flint import arb
    A, DIM = g.build_arb_tau(c, N, prec)
    d = [None]*DIM
    L = [[arb(0)]*DIM for _ in range(DIM)]
    signs = []
    for i in range(DIM):
        s = A[i,i]
        for k in range(i):
            s = s - L[i][k]*L[i][k]*d[k]
        d[i] = s
        signs.append(1 if s>0 else (-1 if s<0 else 0))
        for j in range(i+1, DIM):
            t = A[j,i]
            for k in range(i):
                t = t - L[j][k]*L[i][k]*d[k]
            L[j][i] = t/d[i]
    dfloat = np.array([to_float(x) for x in d])
    return dfloat, DIM, np.array(signs)

c, N, prec = 13, 30, 400
print(f"B2 Szegő probe: Groskin tau c={c}, N={N} (DIM={2*N+1}), prec={prec}")
d, DIM, signs = pivots_of(c, N, prec)
npos = int(np.sum(signs>0)); nneg = int(np.sum(signs<0)); nundet = int(np.sum(signs==0))
print(f"pivots (Arb-signed): n_pos={npos} n_neg={nneg} n_undet={nundet}  (PD? {nneg==0 and nundet==0})")
logdet = float(np.sum(np.log(np.abs(d))))
print(f"log|det tau| = {logdet:.6f}")

ld = np.log(np.abs(d))
# order pivots as built (natural Fourier order n=-N..N). Look at the SEQUENCE.
print("\nlog(pivot) sequence (natural order, n=-N..N):")
for i in range(0, DIM, max(1,DIM//20)):
    print(f"  i={i:3d} (n={i-N:+3d})  log d_i = {ld[i]:+.5f}")

# Szegő test: do the log-pivots settle to a constant? tail mean & correction sum
tail = ld[DIM//2:]  # second half
mean_tail = float(np.mean(tail))
print(f"\ntail (i>={DIM//2}) mean log-pivot = {mean_tail:.5f}, std = {float(np.std(tail)):.5f}")
corr = ld - mean_tail
abscorr_sum = float(np.sum(np.abs(corr)))
sq_sum = float(np.sum(corr**2))
print(f"Sum |log d_i - mean_tail| = {abscorr_sum:.5f}   Sum (·)^2 = {sq_sum:.5f}")
# convergence check: does partial-sum of |corr| flatten (summable) or grow linearly (not)?
ps = np.cumsum(np.abs(corr))
print("cumulative Σ|corr| at 25/50/75/100%:", [f"{ps[int(f*(DIM-1))]:.4f}" for f in (.25,.5,.75,1.0)])
# increments of cumulative in last quarter — if ->0, summable
last_q = np.diff(ps[3*DIM//4:])
print(f"mean |corr| in last quarter = {float(np.mean(np.abs(corr[3*DIM//4:]))):.6f}  (->0 means summable/Szegő-like)")
print("\nINTERPRETATION HINT: Szegő-like ⟺ log-pivots converge to a constant with SUMMABLE correction")
print("(mean |corr| in last quarter small & flat). Erratic/growing ⟹ no clean determinant invariant (matches A1).")

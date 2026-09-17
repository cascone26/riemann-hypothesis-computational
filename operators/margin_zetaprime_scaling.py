#!/usr/bin/env python3
"""
cycle 344 — NEW structural probe: do the CCM/Fejér Weil form's SMALL eigenvalues scale like |zeta'(rho)|^2 ?

§480 derived the NB Gram symbol S(t)=|zeta(1/2+it)|^2/(1/4+t^2) >= 0, vanishing (double) at each zero gamma:
near gamma, |zeta(1/2+it)|^2 ~ |zeta'(1/2+i gamma)|^2 (t-gamma)^2, so S(t) ~ |zeta'(rho)|^2 (t-gamma)^2/(1/4+gamma^2).
PREDICTION: the near-null (small) eigenvalues of the finite Weil form, one per in-window zero, scale like
    mu(gamma) ~ C * |zeta'(1/2+i gamma)|^2 / (1/4+gamma^2).
If true, the operator's small spectrum ENCODES |zeta'| at the zeros => ties positivity margin to zero SIMPLICITY
(|zeta'|!=0) and the Alpoge-Furman 2/3-simple program. REPORT-tier structural test (numerics).
"""
import numpy as np, sys, os, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import l2_fejer as lf
from scipy.linalg import eigh

a, b = 14.0, 44.0          # window with ~8 zeros
N = 20
CEN = np.linspace(a, b, N)
d = lf.build(CEN)
full, G = d['full'], d['overlap']
w, V = eigh(full, G)       # ascending eigenvalues (generalized)
print(f"window [{a},{b}], N={N}: smallest 10 eigenvalues:")
print("  " + "  ".join(f"{x:.3e}" for x in w[:10]))

# in-window zeros
G2000 = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros2000.npy')).astype(float)
zin = G2000[(G2000 >= a) & (G2000 <= b)]
print(f"\nin-window zeros ({len(zin)}): {np.round(zin,3)}")

# |zeta'(1/2+i gamma)|^2 / (1/4+gamma^2)  via mpmath
mp.mp.dps = 25
pred = []
for g in zin:
    zp = mp.zeta(mp.mpc(0.5, g), derivative=1)
    pred.append(float(abs(zp)**2 / (0.25 + g**2)))
pred = np.array(pred)
# match: the smallest len(zin) eigenvalues (sorted asc) vs pred sorted by the SAME ordering?
# both are per-zero; sort both ascending and compare the scaling (ratio should be ~const if law holds)
k = len(zin)
small = np.sort(np.abs(w[:k]))
preds = np.sort(pred)
print(f"\n{'zero gamma':>12} {'|zeta_prime|^2/(1/4+g^2)':>26}")
for g, p in sorted(zip(zin, pred)):
    print(f"{g:>12.3f} {p:>26.5e}")
print(f"\nsmallest {k} |eigenvalues| (asc): {np.array2string(small, formatter={'float':lambda x:f'{x:.3e}'})}")
print(f"pred (asc):                    {np.array2string(preds, formatter={'float':lambda x:f'{x:.3e}'})}")
ratio = small/preds
print(f"ratio eig/pred (asc):          {np.array2string(ratio, formatter={'float':lambda x:f'{x:.3f}'})}")
cv = np.std(ratio)/np.mean(ratio) if np.mean(ratio)!=0 else float('inf')
r_corr = np.corrcoef(np.log(small+1e-300), np.log(preds+1e-300))[0,1]
print(f"\nratio CV = {cv:.3f} (small => law holds, const C); log-log corr = {r_corr:.3f}")
print("VERDICT: if CV small & corr~1, the small eigenvalues scale like |zeta'(rho)|^2/(1/4+g^2) =>")
print("operator small-spectrum encodes |zeta'| at zeros (ties margin to simplicity). REPORT-tier.")

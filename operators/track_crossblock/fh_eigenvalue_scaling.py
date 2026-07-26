#!/usr/bin/env python3
"""B2 / §211c: the collinearity-free Fisher-Hartwig observable — small-eigenvalue index scaling.

Build the certified Groskin tau at a fixed N, get its eigenvalues at HIGH PRECISION (mpmath eigsy;
the smallest are ~e^-35, far below float64 eps ~2e-16, so float64 eigensolve is garbage there).
Test how the SMALL eigenvalues scale with their index k:
  - Fisher-Hartwig symbol zero |theta|^{2a} at DC  ==>  small eigenvalues follow a POWER law  lambda_k ~ C k^{-p}
    (equivalently log lambda_k linear in log k), with p tied to the FH exponent a.
  - A hard gap / exponential tail (lambda_k ~ e^{-c k}) would be a DIFFERENT (non-FH) mechanism.
This is measured at ONE N, so there is no N-vs-logN collinearity (the §211b defeat).

REPORT tier (high-precision numeric).
"""
import importlib.util, os
import numpy as np
import mpmath as mp
from flint import arb

GPATH = os.path.expanduser('~/projects/riemann/external/groskin_2607.02828/code/arb_ldlt_certify.py')
spec = importlib.util.spec_from_file_location('groskin_cert', GPATH)
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)

def to_float(x):
    try: return float(x)
    except Exception:
        try: return float(x.mid())
        except Exception: return float(str(x).split('+/-')[0].strip('[] '))

mp.mp.dps = 220
c, N, prec = 13, 18, 1400
A, DIM = g.build_arb_tau(c, N, prec)
print(f"B2 eigenvalue-scaling: Groskin tau c={c}, N={N}, DIM={DIM}, arb prec={prec}, mp dps={mp.mp.dps}")
# arb midpoints -> mpmath matrix (need the full-precision string, not float64)
def arb_to_mpf(x):
    s = str(x)
    if '+/-' in s: s = s.split('+/-')[0].strip('[] ')
    try: return mp.mpf(s)
    except Exception: return mp.mpf(to_float(x))
Amp = mp.matrix(DIM, DIM)
for i in range(DIM):
    for j in range(DIM):
        Amp[i,j] = arb_to_mpf(A[i,j])

E = mp.eigsy(Amp, eigvals_only=True)
ev = sorted([mp.mpf(e) for e in E])
print(f"min eig = {mp.nstr(ev[0],6)}   max eig = {mp.nstr(ev[-1],6)}   all>0? {all(e>0 for e in ev)}")

# smallest ~half the eigenvalues (the near-null tail), log-log scale
k = np.arange(1, DIM+1)
logl = np.array([float(mp.log(abs(e))) for e in ev])   # sorted ascending
print("\nsmallest eigenvalues (k = 1..):")
for i in range(min(16, DIM)):
    print(f"  k={i+1:2d}  lambda={mp.nstr(ev[i],5):>14}  log(lambda)={logl[i]:+.4f}")

# FH power-law test on the small tail: log lambda_k vs log k  (linear => power law)
m = DIM//2  # use smallest half
lk = np.log(k[:m]); ll = logl[:m]
(p_slope, p_int), *_ = np.linalg.lstsq(np.vstack([lk, np.ones_like(lk)]).T, ll, rcond=None)
pred = p_slope*lk + p_int; rms_pow = float(np.sqrt(np.mean((ll-pred)**2)))
# exponential-tail test: log lambda_k vs k (linear => exponential)
(e_slope, e_int), *_ = np.linalg.lstsq(np.vstack([k[:m], np.ones_like(k[:m])]).T, ll, rcond=None)
predе = e_slope*k[:m] + e_int; rms_exp = float(np.sqrt(np.mean((ll-predе)**2)))

print("\n=== SMALL-EIGENVALUE TAIL (smallest half) ===")
print(f"  POWER law  log l ~ {p_slope:.4f}*log k + {p_int:.4f}   RMS={rms_pow:.4f}   (=> l_k ~ k^{p_slope:.3f})")
print(f"  EXP   tail log l ~ {e_slope:.4f}*k     + {e_int:.4f}   RMS={rms_exp:.4f}   (=> l_k ~ e^({e_slope:.3f} k))")
better = 'POWER LAW (Fisher-Hartwig-like)' if rms_pow < rms_exp else 'EXPONENTIAL (hard-gap-like, NOT FH)'
print(f"  BETTER FIT: {better}   (power RMS {rms_pow:.4f} vs exp RMS {rms_exp:.4f})")
if rms_pow < rms_exp:
    print(f"  => FH-type power law; slope p={p_slope:.3f}. For symbol |theta|^(2a): eigenvalue index-power ~ relates to a.")

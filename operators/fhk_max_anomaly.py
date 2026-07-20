#!/usr/bin/env python3
"""
fhk_max_anomaly.py -- ANOMALY HUNT #2: does max|zeta| in short windows follow the
Fyodorov-Hiary-Keating log-correlated prediction, or deviate?

FHK conjecture: for a window [T, T+2pi],
    max_{t in window} log|zeta(1/2+it)|  ~  loglog T  -  (3/4) logloglog T  +  M
The LEADING coefficient of loglog T is 1; the (3/4) on logloglog T is the sharp
refinement from LOG-CORRELATED-FIELD / branching-random-walk structure (NOT from GUE
pair correlation alone). M is an O(1) fluctuation with a predicted (Gumbel-convolution)
right tail.

Anomaly-capable predictions the data COULD violate:
  (1) leading growth coefficient of E[max] vs loglog T  (predict ~1).
  (2) the logloglog refinement coefficient (predict ~ -0.75).
  (3) the fluctuation distribution shape / right tail (log-correlated => heavier than
      an iid-Gumbel of the same variance; a DIFFERENT tail = a real surprise about the
      organizing structure).
Uses the fast Riemann-Siegel Z(t) evaluator; |zeta(1/2+it)| = |Z(t)|.
"""
import numpy as np
import json, sys, time
sys.path.insert(0, "/Users/scones/projects/riemann/operators")
from deriv_moments_rs import rs_Z_chunk   # vectorized Riemann-Siegel Z(t)

TWO_PI = 2*np.pi
NWIN = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
T_LO, T_HI = 2.0e4, 1.10e6
GRID = 0.03   # sampling step inside each 2pi window (mean spacing ~0.5 at T~1e5)

print(f"[fhk] {NWIN} windows of length 2pi across T in [{T_LO:.0e},{T_HI:.0e}], grid={GRID}", flush=True)
centers = np.linspace(T_LO, T_HI, NWIN)
npts = int(TWO_PI/GRID)
offsets = np.linspace(0, TWO_PI, npts)

t0 = time.time()
maxlogZ = np.empty(NWIN)
CH = 2000   # windows per batch
for i in range(0, NWIN, CH):
    cc = centers[i:i+CH]
    tt = (cc[:, None] + offsets[None, :]).ravel()          # (batch*npts,)
    Z = rs_Z_chunk(tt).reshape(len(cc), npts)
    maxlogZ[i:i+CH] = np.max(np.log(np.abs(Z) + 1e-300), axis=1)
    if (i//CH) % 5 == 0:
        el = time.time()-t0
        print(f"  [fhk] {i+len(cc)}/{NWIN}  {el:.0f}s", flush=True)
print(f"[fhk] done ({time.time()-t0:.0f}s)", flush=True)

# FHK predictors
L1 = np.log(np.log(centers))          # loglog T
L2 = np.log(np.log(np.log(centers)))  # logloglog T

# bin by T, average the max (kills the O(1) fluctuation)
nbin = 24
bins = np.linspace(T_LO, T_HI, nbin+1)
idx = np.digitize(centers, bins) - 1
bmean, bL1, bL2, bstd = [], [], [], []
for b in range(nbin):
    m = idx == b
    if m.sum() < 20: continue
    bmean.append(maxlogZ[m].mean()); bstd.append(maxlogZ[m].std())
    bL1.append(L1[m].mean()); bL2.append(L2[m].mean())
bmean=np.array(bmean); bL1=np.array(bL1); bL2=np.array(bL2); bstd=np.array(bstd)

print("\n" + "="*70)
print("FHK MAXIMUM  E[max log|Z|] = a*loglogT + b*logloglogT + c   (predict a~1, b~-0.75)")
print("="*70)
# fit E[max] ~ a*L1 + b*L2 + c
A = np.vstack([bL1, bL2, np.ones_like(bL1)]).T
coef, *_ = np.linalg.lstsq(A, bmean, rcond=None)
a, b, c = coef
print(f"  leading  a (loglogT)     = {a:.3f}   (FHK predict ~1.0)")
print(f"  refine   b (logloglogT)  = {b:.3f}   (FHK predict ~-0.75)")
# also leading-only fit (a on L1, intercept), since 3-param over narrow lever is noisy
A1 = np.vstack([bL1, np.ones_like(bL1)]).T
a1, c1 = np.linalg.lstsq(A1, bmean, rcond=None)[0]
print(f"  leading-only fit: a={a1:.3f}  (loglogT coefficient, most robust over narrow T-lever)")

# fluctuation distribution: centered max M = maxlogZ - (loglogT - 0.75 logloglogT)
M = maxlogZ - (L1 - 0.75*L2)
mM, sM = M.mean(), M.std()
# right-tail: compare tail decay to a Gumbel of same scale (log-correlated => heavier)
q = np.quantile(M, [0.5, 0.9, 0.99, 0.999])
print(f"\n  centered fluctuation M: mean={mM:.3f} std={sM:.3f}")
print(f"  quantiles [.5,.9,.99,.999] = {np.round(q,3).tolist()}")
print(f"  right-tail width (q99-q50)/std = {(q[2]-q[0])/sM:.3f}  (iid-Gumbel ~2.5; log-corr heavier)")

out = {"n_windows": NWIN, "T_range": [T_LO, T_HI],
       "fit_a_loglogT": float(a), "fit_b_logloglogT": float(b), "fit_c": float(c),
       "leading_only_a": float(a1),
       "fluct_mean": float(mM), "fluct_std": float(sM),
       "fluct_quantiles_50_90_99_999": [float(x) for x in q],
       "right_tail_q99_q50_over_std": float((q[2]-q[0])/sM),
       "binned": {"T_mid": ((bins[:-1]+bins[1:])/2)[:len(bmean)].tolist(),
                  "mean_max": bmean.tolist(), "std_max": bstd.tolist()}}
p = "/Users/scones/projects/riemann/results/fhk_max_anomaly.json"
json.dump(out, open(p, "w"), indent=2)
print(f"\nsaved -> {p}")
print("\nVERDICT: a~1 & b~-0.75 => FHK log-correlated structure confirmed. b far from -0.75,")
print("or tail shape different => a REAL surprise about the max/organizing structure.")
print("(Caveat: T<=1.1e6 gives a NARROW loglog-lever; leading a is robust, b is only indicative.)")

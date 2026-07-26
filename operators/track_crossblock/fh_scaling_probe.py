#!/usr/bin/env python3
"""B2 / §211 follow-up: Fisher-Hartwig vs strong-Szego scaling of the Groskin Weil determinant.

For a range of N, build the CERTIFIED cutoff-free Groskin tau (c fixed), keep the LDL^T in Arb
(prec high enough that ALL pivots are rigorously signed — no float64-before-LDL, the §211 bug),
and record:
  - log|det tau_N| = Σ log d_i     (all d_i > 0 since tau is PD)
  - the well-bottom pivot  min_i log d_i   (the DC near-null: how close to singular)

Strong-Szego:  log det = a*N + c + o(1)           (NO log N term)
Fisher-Hartwig (symbol zero at DC): log det = a*N + b*log N + c,  b = FH exponent (encodes the zero).
Fit both; report b and whether it is nonzero/recognizable. Also fit well-bottom vs N.

REPORT tier (float logs of Arb-signed pivots; PD/inertia itself is Arb-rigorous).
"""
import importlib.util, os
import numpy as np
from flint import arb

GPATH = os.path.expanduser('~/projects/riemann/external/groskin_2607.02828/code/arb_ldlt_certify.py')
spec = importlib.util.spec_from_file_location('groskin_cert', GPATH)
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)

def to_float(x):
    try: return float(x)
    except Exception:
        try: return float(x.mid())
        except Exception: return float(str(x).split('+/-')[0].strip('[] '))

def logdet_and_well(c, N, prec):
    A, DIM = g.build_arb_tau(c, N, prec)
    d = [None]*DIM; L = [[arb(0)]*DIM for _ in range(DIM)]
    undet = 0
    for i in range(DIM):
        s = A[i,i]
        for k in range(i):
            s = s - L[i][k]*L[i][k]*d[k]
        d[i] = s
        if not (s > 0 or s < 0): undet += 1
        for j in range(i+1, DIM):
            t = A[j,i]
            for k in range(i):
                t = t - L[j][k]*L[i][k]*d[k]
            L[j][i] = t/d[i]
    logs = np.array([np.log(abs(to_float(x))) for x in d])
    return float(np.sum(logs)), float(np.min(logs)), undet, DIM

c = 13
prec = 3200
Ns = [8, 10, 12, 14, 16, 18, 20, 22]
rows = []
print(f"FH scaling probe: Groskin tau c={c}, prec={prec}")
for N in Ns:
    ld, well, undet, DIM = logdet_and_well(c, N, prec)
    rows.append((N, DIM, ld, well, undet))
    print(f"  N={N:2d} DIM={DIM:3d}  log|det|={ld:12.5f}  well_bottom={well:10.5f}  undet={undet}", flush=True)

Narr = np.array([r[0] for r in rows], float)
LD = np.array([r[2] for r in rows])
WELL = np.array([r[3] for r in rows])
anyundet = any(r[4] for r in rows)
print(f"\nAll pivots rigorously signed? {not anyundet}  (if False, raise prec)")

# Fit A: strong-Szego  log det = a*N + c
M1 = np.vstack([Narr, np.ones_like(Narr)]).T
(a1, c1), res1, *_ = np.linalg.lstsq(M1, LD, rcond=None)
pred1 = M1 @ [a1, c1]; rms1 = float(np.sqrt(np.mean((LD-pred1)**2)))
# Fit B: Fisher-Hartwig  log det = a*N + b*log N + c
M2 = np.vstack([Narr, np.log(Narr), np.ones_like(Narr)]).T
(a2, b2, c2), *_ = np.linalg.lstsq(M2, LD, rcond=None)
pred2 = M2 @ [a2, b2, c2]; rms2 = float(np.sqrt(np.mean((LD-pred2)**2)))

print("\n=== FIT A (strong-Szego): log det = a*N + c ===")
print(f"  a={a1:.5f}  c={c1:.5f}   RMS residual = {rms1:.5f}")
print("=== FIT B (Fisher-Hartwig): log det = a*N + b*log N + c ===")
print(f"  a={a2:.5f}  b={b2:.5f}  c={c2:.5f}   RMS residual = {rms2:.5f}")
print(f"\n  --> FH log-N coefficient b = {b2:.5f}")
print(f"      Szego RMS {rms1:.4f} vs FH RMS {rms2:.4f}  (FH much better ⟺ real log-N term ⟺ symbol zero)")
print(f"      candidate FH exponents: if b≈-alpha^2, alpha≈{np.sqrt(abs(b2)):.4f} (if b<0)")

# Well-bottom scaling: how fast does the DC near-null approach singular?
Mw = np.vstack([Narr, np.ones_like(Narr)]).T
(wa, wc), *_ = np.linalg.lstsq(Mw, WELL, rcond=None)
Mwl = np.vstack([np.log(Narr), np.ones_like(Narr)]).T
(wla, wlc), *_ = np.linalg.lstsq(Mwl, WELL, rcond=None)
predw = Mw@[wa,wc]; predwl = Mwl@[wla,wlc]
print("\n=== WELL-BOTTOM (DC near-null pivot) scaling ===")
print(f"  linear:  well ≈ {wa:.4f}*N + {wc:.4f}   RMS={float(np.sqrt(np.mean((WELL-predw)**2))):.4f}")
print(f"  log:     well ≈ {wla:.4f}*log N + {wlc:.4f}   RMS={float(np.sqrt(np.mean((WELL-predwl)**2))):.4f}")
print(f"  (linear-in-N deepening ⟹ pivot ~ e^(-const*N), exponential approach to singular = strong FH zero)")

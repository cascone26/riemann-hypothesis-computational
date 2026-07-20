#!/usr/bin/env python3
"""
lehmer_gap_anomaly.py -- ANOMALY HUNT (not a confirmation): does the small-gap tail
of the zeta zeros EXCEED the GUE prediction?

GUE nearest-neighbor spacing (unit mean): p(s) -> (pi^2/3) s^2 as s->0 (quadratic
level repulsion). Sharp, falsifiable predictions the data COULD violate:
  1. small-gap exponent: p(s) ~ C s^alpha, GUE => alpha = 2 exactly.
  2. small-gap coefficient: C = pi^2/3 = 3.2899 (leading), + O(s^4).
  3. cumulative count below eps: N(s<eps)/N ~ (pi^2/9) eps^3  [integral of GUE p(s)].
  4. minimum gap over N zeros: GUE extreme-value => s_min ~ N^{-1/3} (up to log).

A measured EXCESS of small gaps (alpha<2, C too big, or s_min shrinking faster than
N^{-1/3}) = Lehmer-pair excess = zeros coming anomalously close to colliding = a
signature that RH is "barely" true (ties to de Bruijn-Newman Lambda>=0). That would
be a REAL anomaly, not a confirmation. This script measures all four on 2M+ zeros.

Uses raw gamma differences only (no zeta'), so it's fast and exact to the data's 9dp.
"""
import numpy as np
import json, sys

NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 2000000
ZFILE = "/Users/scones/projects/riemann/data/zeros_odlyzko_2M"
TWO_PI = 2*np.pi

print(f"[lehmer] loading {NZ} zeros ...", flush=True)
g = np.loadtxt(ZFILE, max_rows=NZ)
print(f"[lehmer] {len(g)} zeros, range {g[0]:.3f}..{g[-1]:.1f}", flush=True)

# normalized spacings: s_n = (g_{n+1}-g_n) * local_density,
# local mean density of zeros at height t = (1/2pi) log(t/2pi)
gaps = np.diff(g)
mid = 0.5*(g[1:] + g[:-1])
dens = np.log(mid/TWO_PI)/TWO_PI
s = gaps * dens                      # unit-mean normalized spacings
s = s[np.isfinite(s) & (s > 0)]
print(f"[lehmer] {len(s)} normalized spacings, mean={s.mean():.4f} (should be ~1.0)", flush=True)

out = {"n_zeros": len(g), "n_spacings": int(len(s)), "mean_spacing": float(s.mean())}

# ---- (1)&(2) small-gap exponent + coefficient via histogram fit on s in (0, 0.3] ----
edges = np.linspace(0, 0.35, 36)          # fine bins near 0
centers = 0.5*(edges[1:]+edges[:-1])
hist, _ = np.histogram(s, bins=edges, density=True)
fitmask = (centers > 0.03) & (centers < 0.30) & (hist > 0)
# log p = log C + alpha log s
A = np.vstack([np.log(centers[fitmask]), np.ones(fitmask.sum())]).T
alpha, logC = np.linalg.lstsq(A, np.log(hist[fitmask]), rcond=None)[0]
C = np.exp(logC)
print("\n" + "="*70)
print("SMALL-GAP TAIL  p(s) ~ C s^alpha   (GUE: alpha=2, C=pi^2/3=3.290)")
print("="*70)
print(f"  measured alpha = {alpha:.3f}   (GUE 2.000)   {'EXCESS(alpha<2)' if alpha<1.9 else 'consistent' if alpha<2.15 else 'DEFICIT'}")
print(f"  measured C     = {C:.3f}   (GUE 3.290)")
out["small_gap_fit"] = {"alpha": float(alpha), "alpha_gue": 2.0, "C": float(C), "C_gue": float(np.pi**2/3)}

# ---- (3) cumulative small-gap count vs GUE (pi^2/9) eps^3 ----
print("\nCUMULATIVE  frac(s<eps)  vs GUE (pi^2/9) eps^3:")
cum = {}
for eps in [0.05, 0.1, 0.15, 0.2, 0.3]:
    frac = float(np.mean(s < eps))
    gue = (np.pi**2/9)*eps**3
    ratio = frac/gue
    flag = "  <-- EXCESS" if ratio > 1.15 else ("  <-- deficit" if ratio < 0.85 else "")
    print(f"  eps={eps:.2f}: data={frac:.3e}  GUE={gue:.3e}  ratio={ratio:.3f}{flag}")
    cum[str(eps)] = {"data": frac, "gue": gue, "ratio": float(ratio)}
out["cumulative_small_gap"] = cum

# ---- (4) minimum-gap scaling vs GUE N^{-1/3} ----
print("\nMINIMUM normalized gap over first M zeros  (GUE extreme-value: ~ M^{-1/3}):")
mins = {}
for M in [10000, 50000, 200000, 1000000, len(s)]:
    if M > len(s): continue
    smin = float(np.min(s[:M]))
    pred = M**(-1.0/3.0)
    print(f"  M={M:>8}: s_min={smin:.3e}   M^(-1/3)={pred:.3e}   ratio={smin/pred:.3f}")
    mins[str(M)] = {"s_min": smin, "M_pow_-1/3": float(pred), "ratio": float(smin/pred)}
out["min_gap_scaling"] = mins

# smallest few gaps (Lehmer candidates) with their heights
order = np.argsort(s)[:10]
lehmer = [{"rank": int(i+1), "s": float(s[idx]), "gamma": float(g[idx])} for i, idx in enumerate(order)]
print("\n10 SMALLEST normalized gaps (Lehmer candidates):")
for L in lehmer:
    print(f"  s={L['s']:.4e}  at gamma~{L['gamma']:.1f}")
out["lehmer_candidates"] = lehmer

p = "/Users/scones/projects/riemann/results/lehmer_gap_anomaly.json"
json.dump(out, open(p, "w"), indent=2)
print(f"\nsaved -> {p}")
print("\nVERDICT: alpha~2 & ratios~1 => GUE holds (confirmation). alpha<2 or ratios>>1 or")
print("s_min shrinking faster than M^{-1/3} => REAL small-gap EXCESS (Lehmer anomaly).")

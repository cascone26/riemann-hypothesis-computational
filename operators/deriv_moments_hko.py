#!/usr/bin/env python3
"""
deriv_moments_hko.py -- Hughes-Keating-O'Connell / Hughes-Pearce-Crump discrete
derivative-moment test on real verified zeros.

RMT (CUE) PREDICTION (Hughes-Keating-O'Connell 2000; Gonek/Hejhal):
    M_k(T) := (1/N(T)) * sum_{0<gamma_n<=T} |zeta'(1/2 + i gamma_n)|^{2k}
            ~  C_k * (log(T/2pi))^{k(k+2)}
The LEADING POWER k(k+2) is the sharp falsifiable prediction:
    k=1 -> 3,  k=2 -> 8,  k=3 -> 15.
Test: compute M_k at a ladder of T (cumulative over first m zeros), regress
log M_k against L=log(gamma_m/2pi); slope should approach k(k+2).

A DEVIATION from k(k+2) on verified on-line zeros would be genuinely significant.
A match is the strongest fine-grained RMT confirmation on real data (beyond pair
correlation). Uses 31-dp zeros (data/zeros_100k.txt) + mpmath zeta'.
"""
import mpmath as mp
import numpy as np
import json, time, sys

mp.mp.dps = 30
NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
KS = [1, 2, 3]
ZFILE = "/Users/scones/projects/riemann/data/zeros_100k.txt"

print(f"[hko] loading first {NZ} zeros (31 dp) ...", flush=True)
gammas = []
with open(ZFILE) as f:
    for i, line in enumerate(f):
        if i >= NZ:
            break
        gammas.append(mp.mpf(line.strip()))
print(f"[hko] {len(gammas)} zeros, gamma range {float(gammas[0]):.3f}..{float(gammas[-1]):.1f}", flush=True)

# compute |zeta'(1/2+i gamma)|^2 at each zero
absZp2 = np.zeros(len(gammas))
t0 = time.time()
for n, g in enumerate(gammas):
    s = mp.mpf('0.5') + 1j * g
    zp = mp.zeta(s, derivative=1)
    absZp2[n] = float(mp.fabs(zp) ** 2)
    if (n + 1) % 1000 == 0:
        el = time.time() - t0
        print(f"  [hko] {n+1}/{len(gammas)}  {el:.0f}s  eta={el/(n+1)*(len(gammas)-n-1):.0f}s", flush=True)
print(f"[hko] zeta' done ({time.time()-t0:.0f}s)", flush=True)

g = np.array([float(x) for x in gammas])
# cumulative moment ladder
ladder = [500, 1000, 2000, 4000, 6000, 8000, 10000, 15000, 20000, 30000, 50000, len(g)]
ladder = [m for m in ladder if m <= len(g)]
out = {"n_zeros": len(g), "ks": KS, "prediction_power": {k: k*(k+2) for k in KS}, "ladder": {}}
print("\n" + "="*70)
print("HKO DERIVATIVE-MOMENT TEST  (power k(k+2): k=1->3, k=2->8, k=3->15)")
print("="*70)
for k in KS:
    w = absZp2 ** k                      # |zeta'|^{2k}
    Ls, logMs = [], []
    for m in ladder:
        Mk = float(np.mean(w[:m]))
        L = np.log(g[m-1] / (2*np.pi))
        Ls.append(np.log(L)); logMs.append(np.log(Mk))
    Ls = np.array(Ls); logMs = np.array(logMs)
    # slope of log M_k vs log L  == exponent on (log(T/2pi))
    A = np.vstack([Ls, np.ones_like(Ls)]).T
    slope, intercept = np.linalg.lstsq(A, logMs, rcond=None)[0]
    # also local slope over the top half (larger T, cleaner asymptotics)
    half = len(Ls)//2
    sl2, _ = np.linalg.lstsq(A[half:], logMs[half:], rcond=None)[0], None
    slope_top = np.linalg.lstsq(A[half:], logMs[half:], rcond=None)[0][0]
    pred = k*(k+2)
    print(f"\n k={k}: predicted exponent {pred}")
    print(f"   full-range fit slope = {slope:.3f}   top-half slope = {slope_top:.3f}")
    print(f"   M_{k}(T) at largest T = {np.exp(logMs[-1]):.4e}  (L={np.exp(Ls[-1]):.3f})")
    out["ladder"][k] = {"fit_slope_full": float(slope), "fit_slope_top": float(slope_top),
                        "predicted": pred,
                        "M_at_maxT": float(np.exp(logMs[-1])),
                        "logL": Ls.tolist(), "logM": logMs.tolist()}

p = "/Users/scones/projects/riemann/results/deriv_moments_hko.json"
json.dump(out, open(p, "w"), indent=2)
print(f"\nsaved -> {p}")
print("\nNOTE: log-lever is short below T~10^4 (L~0.8..7.4). Slope should TREND toward")
print("k(k+2) from below; clean resolution of exponent 8,15 needs T~10^6+ (millions of zeros).")

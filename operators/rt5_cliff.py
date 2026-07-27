#!/usr/bin/env python3
"""§373: pin the CLIFF — does the Weil-form collapse fill coincide with the exact RvM count?

§372b: the Weil margin cliffs (W/overlap -> 0) at critical density while the overlap frame stays
positive. Locate the cliff precisely: for each band, increase probe count n one at a time; find
n* = the count at which W/overlap first drops below 0.02 (arithmetic collapse turns on). Compare
n* to the exact Riemann-von Mangoldt count of the effective band N(hi)-N(lo), and to the raw zero
count from the actual zero list. If n* == RvM count (per band), §335's "Nyquist in the zeros" is
PINNED to the exact count — a precise, zero-free-computable critical-density law.
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def N_rvm(T):
    if T <= 2 * np.pi: return 0.0
    return T / (2 * np.pi) * np.log(T / (2 * np.pi)) - T / (2 * np.pi) + 7.0 / 8.0

w = 2.0
BANDS = [(14.0, 44.0), (30.0, 60.0), (50.0, 80.0), (75.0, 105.0), (100.0, 130.0),
         (135.0, 165.0), (185.0, 215.0), (245.0, 275.0)]
print("cliff location n* (W/overlap first < 0.02) vs exact RvM count of the effective band [lo-2w, hi+2w]")
print(f"{'band':>10} {'n*':>3} {'RvM_eff':>8} {'rawzeros':>8} {'n*/RvM':>7} {'W@n*-1':>9} {'W@n*':>9}")
for band in BANDS:
    b0, b1 = band
    set_grid(b1 + 30.0)
    lo, hi = b0 - 2 * w, b1 + 2 * w
    rvm_eff = N_rvm(hi) - N_rvm(lo)
    rawz = int(np.sum((L.ZEROS >= lo) & (L.ZEROS <= hi)))
    prev_ratio = None; nstar = None; wprev = None; wat = None
    for n in range(4, 40):
        s = (b1 - b0) / (n - 1)
        cs = np.arange(b0, b1 + 1e-9, s)
        if len(cs) != n: cs = np.linspace(b0, b1, n)
        f = L.build_forms(cs, w)
        mo = float(np.linalg.eigvalsh(f['overlap'])[0])
        mf = float(np.linalg.eigvalsh(f['full'])[0])
        r = mf / mo if mo > 1e-300 else 0.0
        if r < 0.02 and nstar is None:
            nstar = n; wat = mf; wprev = prev_wf
            break
        prev_ratio = r; prev_wf = mf
    if nstar is None:
        print(f"{str([int(b0),int(b1)]):>10}  no cliff < n=40"); continue
    print(f"{str([int(b0),int(b1)]):>10} {nstar:3d} {rvm_eff:8.2f} {rawz:8d} {nstar/rvm_eff:7.3f} "
          f"{wprev:+9.2e} {wat:+9.2e}", flush=True)
print("\nreading: if n*/RvM_eff ~ const (ideally ~1), the arithmetic-collapse fill IS the Weyl count —")
print("§335 pinned: the Weil form cliffs exactly when probe count reaches the zero count of the band,")
print("a precise zero-free critical-density law. Scatter in the constant = the O(1) resolution factor.")

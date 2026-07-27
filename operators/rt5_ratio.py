#!/usr/bin/env python3
"""§372b: does Weil/overlap ratio DEGRADE as a->1? (it MUST — overlap Gram is always PSD, so a
constant ratio to criticality would trivially prove RH). Track the ratio at fills approaching 1,
and — the real point — track it into the regime where the Weil margin can go NEGATIVE while the
overlap floor stays >=0. The ratio's collapse toward/below 0 as a->1 IS the arithmetic content,
localized. Uses the bent-prime control (Lambda(2)x1.05) to show where arithmetic flips sign while
the frame geometry (overlap) cannot.
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def rho(T): return np.log(T / (2 * np.pi)) / (2 * np.pi)

# bent-prime scale (inflate powers of 2 by 1.05) — the t11 sign probe
scale = np.ones_like(L.LAMW)
for i, n in enumerate(L.NS):
    m = n
    while m % 2 == 0: m //= 2
    if m == 1: scale[i] = 1.05

w = 2.0
print("Weil/overlap ratio approaching criticality; bent-prime margin (arith sign flip) alongside.")
print("overlap Gram is ALWAYS PSD (>=0); Weil & bent can go negative — that gap IS the arithmetic.")
for band in [(14.0, 44.0), (100.0, 130.0)]:
    set_grid(band[1] + 30.0)
    Tc = (band[0] + band[1]) / 2
    Nz = rho(Tc) * (band[1] - band[0]) + 1.0
    print(f"\nband {band}  Weyl~{Nz:.1f}")
    print(f"{'a':>5} {'n':>3} {'overlap>=0':>10} {'Weil':>11} {'W/ovl':>7} {'bent':>11} {'bent<0?':>7}")
    for a in (0.7, 0.8, 0.85, 0.9, 0.93, 0.96, 0.99, 1.02, 1.05):
        ntarget = max(4, int(round(a * Nz)))
        s = (band[1] - band[0]) / (ntarget - 1)
        cs = np.arange(band[0], band[1] + 1e-9, s)
        f = L.build_forms(cs, w)
        fb = L.build_forms(cs, w, lam_scale=scale)
        mo = float(np.linalg.eigvalsh(f['overlap'])[0])
        mf = float(np.linalg.eigvalsh(f['full'])[0])
        mb = float(np.linalg.eigvalsh(fb['full'])[0])
        r = mf / mo if mo > 1e-300 else float('nan')
        print(f"{a:5.2f} {len(cs):3d} {mo:10.2e} {mf:+11.2e} {r:7.3f} {mb:+11.2e} "
              f"{'NEG' if mb < 0 else 'pos':>7}", flush=True)
print("\nreading: if W/ovl stays ~0.5 to a=1 and beyond, the arithmetic is a benign bounded")
print("perturbation (the §372 reduction is tight). If W/ovl DROPS toward 0 as a->1 while overlap")
print("stays positive, the arithmetic is precisely the ratio-deficit at criticality — and the")
print("bent-prime going NEG where the true stays >=0 marks the exact sign-sensitive corner (§335).")

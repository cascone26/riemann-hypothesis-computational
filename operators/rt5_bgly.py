#!/usr/bin/env python3
"""§371: does OUR margin plunge match the BGLy (1-a^2)^2 frame-bound rate (§370)?

BGLy 2009: Gaussian Gabor lower frame bound A(a) ~ (1-a^2)^2 as density a -> 1^-.
Our density a = (#probes)/(Weyl count Nz) = probe density / critical density.
For each band, sweep s (=> a from ~0.5 up to 1), fit log(true margin) vs log(1-a^2):
BGLy predicts SLOPE 2. Also compare true-margin vs overlap-Gram lmin (the pure L2/frame
object, no arithmetic) — BGLy is about the FRAME (overlap) bound, so the overlap Gram lmin
should show the cleaner (1-a^2)^2; the full Weil margin may inherit it or deviate (the deviation
= the arithmetic content on top of the frame geometry).
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def rho(T): return np.log(T / (2 * np.pi)) / (2 * np.pi)

print("BGLy check: log(margin) vs log(1-a^2), a = #probes/Weyl-count.  BGLy predicts slope 2.")
print(f"{'band':>10} {'w':>4} | {'slope full':>11} {'slope overlap':>14} (BGLy=2)")
for band in [(14.0, 44.0), (50.0, 80.0), (100.0, 130.0)]:
    set_grid(band[1] + 30.0)
    Tc = (band[0] + band[1]) / 2
    Nz = rho(Tc) * (band[1] - band[0]) + 1.0     # approx Weyl count in the window
    for w in (2.0, 3.0):
        A_full, A_ovl, X = [], [], []
        for s in np.arange(6.0, 1.05, -0.05):
            cs = np.arange(band[0], band[1] + 1e-9, s)
            n = len(cs)
            if n < 4: continue
            a = n / Nz
            if a >= 0.999: break
            f = L.build_forms(cs, w)
            mf = float(np.linalg.eigvalsh(f['full'])[0])
            mo = float(np.linalg.eigvalsh(f['overlap'])[0])
            if mf < 1e-13 or mo < 1e-13: break
            if 0.55 < a < 0.98:
                X.append(np.log(1 - a * a)); A_full.append(np.log(mf)); A_ovl.append(np.log(mo))
        if len(X) < 4:
            print(f"{str([int(band[0]),int(band[1])]):>10} {w:4.1f} | too few points"); continue
        X = np.array(X)
        sf = np.polyfit(X, np.array(A_full), 1)[0]
        so = np.polyfit(X, np.array(A_ovl), 1)[0]
        print(f"{str([int(band[0]),int(band[1])]):>10} {w:4.1f} | {sf:11.2f} {so:14.2f}", flush=True)
print("\nreading: overlap-Gram slope ~2 confirms our collapse IS the BGLy frame-bound plunge (pure")
print("geometry). If the FULL Weil margin slope differs from 2, that gap is the arithmetic content")
print("riding on the frame collapse — the quantity the certificate must control beyond frame theory.")

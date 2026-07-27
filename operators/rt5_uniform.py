#!/usr/bin/env python3
"""§372: the UNIFORM SUB-CRITICAL MARGIN — empirical shadow of the open Q1 bound (§370).

The frame literature (BGLy) has no config-independent lower bound on lambda_min below critical
density. Our question: at FIXED sub-critical fill a = #probes/Weyl-count (a in {0.6,0.7,0.8,0.9}),
is the TRUE Weil margin bounded below UNIFORMLY across bands T=29..440, and how does that floor
depend on w? A positive, roughly-band-independent floor at fixed a<1 is exactly the object a
uniform theorem would assert — measured. Also track the OVERLAP-Gram floor (pure frame) alongside,
to see whether the Weil (arithmetic) margin's uniform floor is better/worse than the frame one.
Fill a is set by choosing s so that #probes = round(a * Weyl-count).
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def rho(T): return np.log(T / (2 * np.pi)) / (2 * np.pi)

CENTERS = [29, 45, 65, 90, 115, 150, 200, 260, 340, 440]
WID = 30.0
FILLS = [0.6, 0.7, 0.8, 0.9]
for w in (3.0, 2.0):
    print(f"\n=== w={w}: TRUE Weil margin at fixed sub-critical fill a (overlap-Gram floor in parens) ===")
    print(f"{'band':>10} {'Nz':>5} | " + " | ".join(f"a={a}" for a in FILLS))
    floors = {a: [] for a in FILLS}
    ofloors = {a: [] for a in FILLS}
    for Tc in CENTERS:
        b0, b1 = Tc - WID / 2, Tc + WID / 2
        set_grid(b1 + 30.0)
        Nz = rho(Tc) * WID + 1.0
        cells = []
        for a in FILLS:
            ntarget = max(4, int(round(a * Nz)))
            s = WID / (ntarget - 1)
            cs = np.arange(b0, b1 + 1e-9, s)
            f = L.build_forms(cs, w)
            mf = float(np.linalg.eigvalsh(f['full'])[0])
            mo = float(np.linalg.eigvalsh(f['overlap'])[0])
            floors[a].append(mf); ofloors[a].append(mo)
            cells.append(f"{mf:8.1e}({mo:7.1e})")
        print(f"{str([int(b0),int(b1)]):>10} {Nz:5.1f} | " + " | ".join(cells), flush=True)
    print(f"  --- w={w} UNIFORM FLOOR (min over bands) ---")
    for a in FILLS:
        mn, mx = min(floors[a]), max(floors[a])
        omn = min(ofloors[a])
        print(f"    a={a}: Weil inf={mn:.2e}  sup={mx:.2e}  (ratio {mx/mn:.1f})  | overlap inf={omn:.2e}",
              flush=True)
print("\nreading: a positive Weil inf-floor at fixed a<1, roughly band-flat (ratio ~O(1)), is the")
print("empirical form of the OPEN uniform sub-critical bound (§370). If the Weil floor is band-flat")
print("where the overlap floor is too, the arithmetic rides the frame geometry uniformly — the")
print("uniform bound would then reduce to the (open but concrete) frame-theoretic one.")

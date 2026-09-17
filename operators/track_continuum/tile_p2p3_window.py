#!/usr/bin/env python3
"""
cycle 342 — continuum SOS tiling of the p={2,3} window [0.549, 0.567] (even sector, modes [1,3,5,7]).

Upgrades cy339-341's discrete point-certs (a=0.549,0.565,0.567) to an INTERVAL claim: does the even-sector
Weil form Q(a) ⪰ 0 for ALL a in each tile, via SOS on the affine interval-enclosure of the matrix entries?
REPORT-TIER (float64 SDP per sos_affine_compose header) — a certified tile = strong evidence of interval
positivity; a fail = the knife-edge blocks the float continuum (then Arb interval-SOS or the analytic relay).
"""
import numpy as np
from tile_sos_continuum import tile_and_certify

# fine tiles across the p={2,3} sub-window (narrow => tighter affine enclosure => best chance at the knife-edge)
edges = np.round(np.arange(0.549, 0.5675, 0.003), 4)
windows = [(float(edges[i]), float(edges[i+1])) for i in range(len(edges)-1)]
print(f"Tiling p={{2,3}} window into {len(windows)} tiles: {windows}")
results = tile_and_certify(windows, verbose=True)
print("\n=== SUMMARY ===")
n_pass = sum(1 for r in results if r.get('passed'))
for r in results:
    print(f"  [{r['a0']:.4f},{r['a1']:.4f}] w={r['width']:.4f}: {r['status']}  gram_min_eig={r.get('gram_min_eig')}")
print(f"\n{n_pass}/{len(results)} tiles certified. "
      + ("Even-sector INTERVAL positivity across the certified tiles (Report-tier)."
         if n_pass else "No tile certified — knife-edge blocks the float continuum; needs Arb interval-SOS / analytic relay."))

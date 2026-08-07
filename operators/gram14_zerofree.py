#!/usr/bin/env python3
"""14x14 zero-free Weil Gram on probes 12..46 step 2, 4W=1.32 (the ambient basis for the top-9
concentration eigenspace pencil cert — anchor assembly, §804). Saves mids+rads to gram14.json."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
fourW = 1.32; RC = 400
CENTERS = [12+2*i for i in range(18) if 12+2*i <= 44]  # 12..44 step 2 -> 17 probes; cap at BAND1=44
ctx.prec = 220
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
cs = [arb(c) for c in CENTERS]; n = len(cs)
print(f"{n}x{n} zero-free Gram on {CENTERS} at 4W={fourW}...", flush=True)
Q = {}
for j in range(n):
    for k in range(j, n):
        e = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(RC)) - ac.prime_entry(cs[j], cs[k])
        Q[(j,k)] = (float(e.mid()), float(e.rad()))
    print(f"  row {j}/{n-1} done", flush=True)
G = {}
for j in range(n):
    for k in range(j, n):
        e = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))
        G[(j,k)] = (float(e.mid()), float(e.rad()))
out = {"centers": CENTERS, "fourW": fourW, "RC": RC,
       "Q": {f"{j},{k}": v for (j,k), v in Q.items()},
       "G": {f"{j},{k}": v for (j,k), v in G.items()}}
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gram14.json"), "w").write(json.dumps(out))
print("saved gram14.json", flush=True)

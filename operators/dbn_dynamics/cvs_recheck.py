#!/usr/bin/env python3
"""§169 recheck: resolve the diag c=21 monotonicity anomaly at higher resolution."""
import sys, json, time
from multiprocessing import Pool
import mpmath as mp
import connes_cvs.operator as O
from cvs_interval import build_decoupled
N, dps, T = 120, 250, 600
GRID = [(19,19),(21,21),(23,23)]  # c=19 reference (expect ~-90), c=21 (expect deepen to ~-95+), c=23
def run(cc):
    cL,cP=cc; t0=time.time(); mp.mp.dps=dps
    Q,L=build_decoupled(cL,cP,N,T,dps)
    lam,vec=O.compute_ground_state(Q)
    r={"c":cL,"N":N,"dps":dps,"T":T,"log10_lam":(float(mp.log10(abs(lam))) if lam!=0 else None),"secs":round(time.time()-t0,1)}
    print(f"c={cL}: log10|lam|={r['log10_lam']}  ({r['secs']}s)",flush=True); return r
if __name__=="__main__":
    print(f"§169 recheck diag c=19,21,23 at N={N} dps={dps} T={T}, flint={O.HAS_FLINT}",flush=True)
    with Pool(3) as p: res=p.map(run,GRID)
    json.dump(res,open("cvs_recheck_result.json","w"),indent=2)
    print("\nRECHECK:"); [print(f"  c={r['c']}: log10|lam|={r['log10_lam']}") for r in sorted(res,key=lambda r:r['c'])]
    print("DONE")

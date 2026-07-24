#!/usr/bin/env python3
"""§170: dense between-prime map on gap [19,23) — validate §169b two-scale structure & fit the deepening law."""
import sys, json, time
from multiprocessing import Pool
import mpmath as mp
import connes_cvs.operator as O
from cvs_interval import build_decoupled
N, dps, T = 124, 260, 600
# gap [19,23): fix complete primes<=19, sweep interval cL; then gate test at 23 (cP=19 missing 23) + restore (cP=23)
GRID = [("gap",19.0,19),("gap",20.0,19),("gap",21.0,19),("gap",22.0,19),("gap",22.9,19),
        ("gate",23.0,19),("restore",23.0,23)]
def run(a):
    tag,cL,cP=a; t0=time.time(); mp.mp.dps=dps
    Q,L=build_decoupled(cL,cP,N,T,dps); lam,_=O.compute_ground_state(Q)
    r={"tag":tag,"cL":cL,"cP":cP,"L":float(L),"log10_lam":(float(mp.log10(abs(lam))) if lam!=0 else None),"secs":round(time.time()-t0,1)}
    print(f"[{tag}] cL={cL} cP={cP} (L={r['L']:.4f}): log10|lam|={r['log10_lam']}  ({r['secs']}s)",flush=True); return r
if __name__=="__main__":
    print(f"§170 between-prime map [19,23) — N={N} dps={dps} T={T}, flint={O.HAS_FLINT}",flush=True)
    with Pool(7) as p: res=p.map(run,GRID)
    json.dump(res,open("cvs_between_result.json","w"),indent=2)
    gap=sorted([r for r in res if r['tag']=='gap'],key=lambda r:r['cL'])
    print("\n[gap] cP=19 fixed, sweep interval (want smooth concave deepening):")
    for i,r in enumerate(gap):
        s=("" if i==0 else f"  slope={ (r['log10_lam']-gap[i-1]['log10_lam'])/(r['L']-gap[i-1]['L']):.1f}/L")
        print(f"   cL={r['cL']} (L={r['L']:.4f}): log10|lam|={r['log10_lam']}{s}")
    gate=[r for r in res if r['tag']=='gate'][0]; rest=[r for r in res if r['tag']=='restore'][0]
    print(f"\n[gate]    cL=23 cP=19 (missing prime 23): log10|lam|={gate['log10_lam']}  (expect COLLAPSE ~O(1))")
    print(f"[restore] cL=23 cP=23 (complete):          log10|lam|={rest['log10_lam']}  (expect near-null restored, deeper than 22.9)")
    print("DONE")

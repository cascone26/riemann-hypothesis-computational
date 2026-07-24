#!/usr/bin/env python3
"""§172 §6.6(a): is the CvS smallest eigenvalue SIMPLE (well-separated) and EVEN as c grows? Tests two of
Connes' three open requirements. Computes lowest-4 of the FULL operator + ground-vector parity, vs c."""
import sys, json, time
from multiprocessing import Pool
import mpmath as mp
import connes_cvs.operator as O
from cvs_interval import build_decoupled
N, dps, T = 90, 180, 400
CS = [13, 15, 17, 19]
def run(c):
    t0=time.time(); mp.mp.dps=dps
    Q,L = build_decoupled(c, c, N, T, dps)     # standard operator (cL=cP=c)
    ev, V = mp.eigsy(Q)                          # FULL operator, both sectors
    order = sorted(range(Q.rows), key=lambda i: ev[i])
    lo = [ev[order[k]] for k in range(4)]
    # parity of ground eigenvector: even => v[N+k] == v[N-k]
    g = [V[r, order[0]] for r in range(Q.rows)]
    asym = sum(abs(g[N+k]-g[N-k]) for k in range(1,N+1)) / (sum(abs(x) for x in g)+mp.mpf(1e-300))
    r={"c":c,"lam0":float(mp.log10(abs(lo[0]))),"lam1":float(mp.log10(abs(lo[1]))),
       "lam2":float(mp.log10(abs(lo[2]))),"gap_ratio_log10":float(mp.log10(abs(lo[1])/abs(lo[0]))),
       "ground_asym":float(asym),"secs":round(time.time()-t0,1)}
    print(f"c={c}: log10|λ0|={r['lam0']:.2f} |λ1|={r['lam1']:.2f} gap(λ1/λ0)=10^{r['gap_ratio_log10']:.1f} "
          f"ground_asymmetry={r['ground_asym']:.2e} ({r['secs']}s)",flush=True)
    return r
if __name__=="__main__":
    print(f"§172 simplicity+parity — N={N} dps={dps} T={T}, flint={O.HAS_FLINT}",flush=True)
    with Pool(4) as p: res=p.map(run,CS)
    json.dump(res,open("cvs_simple_result.json","w"),indent=2)
    print("\nSIMPLE if gap(λ1/λ0) stays large (many OOM); EVEN if ground_asymmetry ~ 0.")
    for r in sorted(res,key=lambda r:r['c']):
        print(f"  c={r['c']}: gap=10^{r['gap_ratio_log10']:.1f}  asym={r['ground_asym']:.1e}")
    print("DONE")

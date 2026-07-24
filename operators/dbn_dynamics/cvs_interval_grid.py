#!/usr/bin/env python3
"""
cvs_interval_grid.py — LEDGER §169 grid driver (Olympus, 16-core, flint). Runs the interval-vs-prime
decoupling grid in parallel and writes JSON. Each (cL,cP) point is independent → multiprocessing.Pool.

Grid:
  (A) FIX primes ≤13, EXTEND interval: cL ∈ {13,14,15,17,19,21}, cP=13.  (isolate interval effect)
  (B) FIX interval log21, ADD primes:  cL=21, cP ∈ {13,17,19,21}.        (isolate prime effect)
Measures λ_min and γ1 error at each point. c-limited N (default 90, dps 200, T 400).

Usage: cvs_interval_grid.py [N] [dps] [T] [nproc]
"""
import sys, json, time
from multiprocessing import Pool
import mpmath as mp
import connes_cvs.operator as O
from cvs_interval import build_decoupled

N     = int(sys.argv[1]) if len(sys.argv) > 1 else 100
dps   = int(sys.argv[2]) if len(sys.argv) > 2 else 200
T     = int(sys.argv[3]) if len(sys.argv) > 3 else 400
nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 11

# Decomposition grid. Primary metric = log10|λ_min| (bit-validated vs Groskin; |γ1 err| ∝ λ_min, ratio ~13000).
#   diag : standard c-sweep cL=cP  (the real convergence sequence, Groskin-comparable)
#   A    : interval-only — FIX primes ≤13, extend interval cL  (isolate the geometric effect)
#   B    : prime-only    — FIX interval log21, add primes cP    (isolate the arithmetic effect)
_pts = {
    ("diag", 13, 13), ("diag", 17, 17), ("diag", 19, 19), ("diag", 21, 21),
    ("A", 14, 13), ("A", 15, 13), ("A", 17, 13), ("A", 19, 13), ("A", 21, 13),
    ("B", 21, 17), ("B", 21, 19),
}
GRID = sorted(_pts, key=lambda x: (x[0], x[1], x[2]))

def run_point(args):
    tag, cL, cP = args
    t0 = time.time()
    mp.mp.dps = dps
    Q, L = build_decoupled(cL, cP, N, T, dps)
    lam, vec = O.compute_ground_state(Q)
    zs = O.extract_zeros(vec, float(L), n_zeros=1, dps=dps)
    e1 = zs[0]['error']
    rec = {
        "tag": tag, "cL": cL, "cP": cP, "N": N, "dps": dps, "T": T,
        "L": float(L),
        "lambda_min": mp.nstr(lam, 12),
        "log10_lam": (float(mp.log10(abs(lam))) if lam != 0 else None),
        "g1_error": (mp.nstr(e1, 8) if e1 is not None else None),
        "log10_g1err": (float(mp.log10(e1)) if (e1 is not None and e1 > 0) else None),
        "secs": round(time.time() - t0, 1),
    }
    print(f"[{tag}] cL={cL} cP={cP}: log10|λ|={rec['log10_lam']}, "
          f"log10 γ1err={rec['log10_g1err']}  ({rec['secs']}s)", flush=True)
    return rec

if __name__ == "__main__":
    print(f"§169 interval/prime decoupling grid — N={N} dps={dps} T={T} nproc={nproc}, "
          f"{len(GRID)} points, flint={O.HAS_FLINT}", flush=True)
    t0 = time.time()
    with Pool(nproc) as pool:
        results = pool.map(run_point, GRID)
    out = {"meta": {"N": N, "dps": dps, "T": T, "elapsed": round(time.time()-t0, 1)}, "points": results}
    with open("cvs_interval_grid_result.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nDONE in {out['meta']['elapsed']}s -> cvs_interval_grid_result.json")
    def row(r): return f"log10|λ|={r['log10_lam']}  (log10 γ1err={r['log10_g1err']})"
    diag = sorted([r for r in results if r["tag"]=="diag"], key=lambda r: r["cL"])
    A    = sorted([r for r in results if r["tag"]=="A"] + [r for r in diag if r["cL"]==13], key=lambda r: r["cL"])
    B    = sorted([r for r in results if r["tag"]=="B"] + [r for r in diag if r["cL"]==21], key=lambda r: r["cP"])
    print("\n[diag] standard c-sweep cL=cP (Groskin-comparable convergence):")
    for r in diag: print(f"   c={r['cL']:>3} (L={r['L']:.4f}): {row(r)}")
    print("\n[A] INTERVAL-ONLY: fix primes≤13, extend interval cL — does log10|λ| keep dropping? (⟹ geometry drives it)")
    for r in A: print(f"   cL={r['cL']:>3} cP=13 (L={r['L']:.4f}): {row(r)}")
    print("\n[B] PRIME-ONLY: fix interval log21, add primes cP — small change? (⟹ arithmetic doesn't drive it)")
    for r in B: print(f"   cL=21 cP={r['cP']:>3}: {row(r)}")
    print("\nDECODE: compare each diag step's drop vs the interval-only(A) drop over the same ΔL.")
    print("If A reproduces most of diag's drop and B is flat ⟹ Groskin (e) confirmed with full control: geometric convergence.")

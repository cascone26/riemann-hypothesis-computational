#!/usr/bin/env python3
"""B2 (now the primary line): does the LOCAL symmetric root √τ converge to a fixed limit operator M∞ as
c→∞? B1b showed √τ is local & c-stable. If its band entries M[n,n+r], viewed as functions of frequency
ω=2πn/L, converge to c-independent profiles, then √τ → a fixed local (pseudodifferential) operator M∞ with
a symbol m(ω,r), and τ∞ = M∞·M∞ ⪰ 0 is MANIFESTLY positive — a proof shape that survives the infinite zero
tower (a bounded local operator acts uniformly at all frequencies). This tests that convergence.

√τ is dominated by O(1) eigenvalues so double-precision eigh on the Arb-assembled τ suffices for band
entries. We compare the r=0,1,2 subdiagonal profiles M[n,n+r](ω) across c on a shared ω-grid.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import numpy as np
from flint import arb, ctx
from arb_ldlt_certify import build_arb_tau

def sqrt_tau(c,N,prec=2000):
    ctx.prec=prec; A,DIM=build_arb_tau(c,N,prec)
    T=np.zeros((DIM,DIM))
    for i in range(DIM):
        for j in range(i,DIM):
            v=float(A[i,j].mid()); T[i,j]=v; T[j,i]=v
    w,V=np.linalg.eigh(T); w=np.clip(w,0,None)
    R=(V*np.sqrt(w)[None,:])@V.T
    return R,DIM

def profiles(R,N,L):
    """M[n,n+r] as (omega, value) for r=0,1,2 (omega=2pi n/L)."""
    DIM=R.shape[0]; out={}
    for r in (0,1,2):
        pts=[]
        for i in range(DIM-r):
            n=i-N; om=2*math.pi*n/L
            pts.append((om, R[i,i+r]))
        out[r]=pts
    return out

def interp(pts, xq):
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    return np.interp(xq, xs, ys)

def main():
    N=80
    cs=[50,200,800,3200]
    print(f"=== B2: convergence of √τ band symbol M[n,n+r](ω) across c (N={N}) ===")
    prof_by_c={}
    for c in cs:
        R,DIM=sqrt_tau(c,N); L=math.log(c)
        prof_by_c[c]=(profiles(R,N,L),L)
        print(f"  built √τ at c={c} (L={L:.3f})")
    # compare on shared omega grid [-8,8] (the interesting band incl. neg region and near-null)
    xq=np.linspace(-8,8,33)
    print("\n  r  omega  " + "  ".join(f"c={c}" for c in cs))
    conv={}
    for r in (0,1,2):
        rows=[]
        for c in cs:
            pr,L=prof_by_c[c]; rows.append(interp(pr[r],xq))
        rows=np.array(rows)
        # convergence: max |c=3200 - c=800| over the grid, relative to typical magnitude
        drift=float(np.max(np.abs(rows[-1]-rows[-2])))
        scale=float(np.mean(np.abs(rows[-1]))) or 1
        conv[r]=drift/scale
        print(f"\n  r={r}: last-step drift/scale = {drift/scale:.3f}")
        for k in range(0,33,4):
            print(f"    ω={xq[k]:+5.1f}  " + "  ".join(f"{rows[ci,k]:+7.4f}" for ci in range(len(cs))))
    json.dump({"omega_grid":xq.tolist(),
               "profiles":{str(c):{str(r):prof_by_c[c][0][r] for r in (0,1,2)} for c in cs},
               "convergence":conv},
              open("results/sqrt_tau_limit_B2.json","w"),indent=2)
    print("\n=== VERDICT (B2) ===")
    worst=max(conv.values())
    print(f"  worst r-profile last-step drift/scale = {worst:.3f}")
    if worst<0.10:
        print("  => √τ band symbol CONVERGES: √τ → a fixed local operator M∞(ω). τ∞=M∞² manifestly PSD.")
        print("     This is a proof-shaped limit that survives the infinite zero tower. Pursue M∞'s")
        print("     closed form (its symbol) → attempt the operator-level positivity proof (rung 4-5).")
    elif worst<0.30:
        print("  => √τ symbol nearly converges (mild drift); larger c / better norm may confirm. Promising.")
    else:
        print("  => √τ symbol does NOT converge cleanly; the local factor drifts with c. Reassess Vector B.")
    print("wrote results/sqrt_tau_limit_B2.json")

if __name__=="__main__":
    main()

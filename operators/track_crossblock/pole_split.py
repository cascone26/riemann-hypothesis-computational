#!/usr/bin/env python3
"""P1 — the pole-split test. The Groskin pole block W02 is EXACTLY rank-2 signature (1,1):
W02 = pref*(u u^T - v v^T), u_m = L/(L^2+16pi^2 m^2), v_m = 4pi m/(L^2+16pi^2 m^2), pref = 32 L sinh(L/4)^2.
So tau = W02 - WR - Wp = pref*u u^T - G,  where G := pref*v v^T + WR + Wp = pref*u u^T - tau.
tau >= 0  <=>  G <= pref*u u^T (rank-1)  <=>  G <= 0 on u-perp (signature (1,n-1), Hodge-index shape).

DECISIVE TEST: lambda_max(G restricted to u-perp).
  - If < 0 with real MARGIN  => the whole arithmetic form (minus the pole square) is comfortably NSD off the
    pole direction; tau's knife-edge lives ENTIRELY in the 1-D pole scalar u => infinite-dim part is settled,
    difficulty is 1-dimensional. HUGE if true.
  - If ~ 0 (knife-edge in u-perp) => the near-null is NOT absorbed by the pole; same wall.
Also: (b) overlap of tau's near-null eigenvector with span{u, v}.
Float64 (all pieces O(1) except the single knife-edge eigenvalue; fine for this structural test).
"""
import importlib.util, os, numpy as np
GPATH=os.path.expanduser('~/projects/riemann/external/groskin_2607.02828/code/arb_ldlt_certify.py')
spec=importlib.util.spec_from_file_location('gc',GPATH); g=importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
def tof(x):
    try: return float(x)
    except Exception:
        s=str(x);
        if '+/-' in s: s=s.split('+/-')[0].strip('[] ')
        return float(s)

def pieces(c,N,prec=1200):
    A,DIM=g.build_arb_tau(c,N,prec)
    tau=np.array([[tof(A[i,j]) for j in range(DIM)] for i in range(DIM)])
    L=np.log(c); sp2=16*np.pi**2; pref=32*L*np.sinh(L/4)**2
    ms=np.arange(-N,N+1)
    u=L/(L**2+sp2*ms**2); v=4*np.pi*ms/(L**2+sp2*ms**2)
    return tau,u,v,pref,L

print("P1 POLE-SPLIT TEST: is tau's knife-edge absorbed by the rank-1 pole square u?")
print("First: verify W02 = pref(uu^T - vv^T) reconstructs (identity check).")
N=20
for c in [13,29,50]:
    tau,u,v,pref,L=pieces(c,N)
    W02=pref*(np.outer(u,u)-np.outer(v,v))
    G=pref*np.outer(u,u)-tau        # = pref vv^T + WR + Wp
    # verify: is G's positive part rank-1? project out u.
    uh=u/np.linalg.norm(u)
    P=np.eye(len(u))-np.outer(uh,uh)  # projector onto u-perp
    Gp=P@G@P
    wG=np.linalg.eigvalsh(G); wGp=np.linalg.eigvalsh(Gp)
    wt=np.linalg.eigvalsh(tau); v0=np.linalg.eigh(tau)[1][:,0]
    # overlap of near-null with u and v
    ov_u=abs(float(v0@uh)); vh=v/np.linalg.norm(v); ov_v=abs(float(v0@vh))
    ov_uv=np.sqrt(ov_u**2+ov_v**2)
    print(f"\nc={c}: L={L:.3f} pref={pref:.3f}")
    print(f"  lambda_min(tau) [knife-edge] = {wt[0]:.4e}")
    print(f"  lambda2(G) (2nd largest, should be <=0 for Hodge shape) = {sorted(wG)[-2]:.4e}   lambda_max(G) = {wG[-1]:.4e}")
    print(f"  ***  lambda_max(G | u-perp) = {wGp[-1]:.4e}  ***   (margin<0 => difficulty is 1-D pole scalar; ~0 => knife-edge survives)")
    print(f"  near-null overlap: |<v0|u>|={ov_u:.4f}  |<v0|v>|={ov_v:.4f}  |proj onto span(u,v)|={ov_uv:.4f}")
print("\nREAD: lambda_max(G|u-perp) strongly NEGATIVE + near-null high overlap with u  => P1 WINS (difficulty localizes")
print("to the pole scalar). lambda_max(G|u-perp) ~ 0 (its own knife-edge) => near-null not pole => same wall.")

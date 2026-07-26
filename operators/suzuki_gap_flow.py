#!/usr/bin/env python3
"""suzuki_gap_flow.py -- §241: the a-FLOW of the Suzuki spectrum: branch trajectories + gap persistence.

NEW FRAME (Fable brainstorm 2026-07-24 late): RH + simple zeros <=> for every fixed k, the k-th spectral
gap of D_{a,theta} stays bounded below as a->infinity. An off-line zero pair can only hide in the (proven
real, Thm 1.5) spectra as TWO BRANCHES COLLIDING IN THE LIMIT at the pair's common ordinate. So the object
to study is the FLOW z_k(a) and gap_k(a) = z_{k+1}(a)-z_k(a):
  - between thresholds a = (1/2)log q (prime powers q): smooth scaling flow;
  - at each threshold: a shift-operator coupling enters CUBICALLY SOFTLY (Onset Lemma §219);
  - prediction to test: branches deform from the regular momentum lattice (spacing pi/a-ish) toward the
    IRREGULAR zeta gaps, driven by the prime insertions; no finite-a collisions (interlacing/simplicity),
    gaps heading to |gamma_{k+1}-gamma_k|.
MEASURE: (1) branch continuity/smoothness across thresholds 0.347(2), 0.549(3), 0.693(4), 0.805(5), 0.973(7);
(2) min over k<=10 of gap_k(a) -- the collision detector; (3) per-branch drift rate d z_k/da between vs at
thresholds (softness check); (4) approach of low gaps to the zeta gaps 6.89, 3.99, 5.41, 2.51, 4.65...
theta protocol: fixed theta in {0, pi/2, pi, 3pi/2} (each a self-adjoint family; branch ID well-defined).
Reuses suzuki_W_zeros machinery verbatim (M=2000,K=50 for speed; docstring warns a>~1 discretization).
"""
import numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_W_zeros as SW

THRESH = {0.3466: 2, 0.5493: 3, 0.6931: 4, 0.8047: 5, 0.9730: 7}
GAMMAS = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862, 40.9187, 43.3271, 48.0052, 49.7738]

def zeros_at(a, theta, M=2000, K=50, zmax=60.0, ng=6000):
    vhp, vhm, lam_a = SW.v_hat_pm(a, M=M, K=K)
    zg = np.linspace(0.5, zmax, ng)
    Wz = np.abs(np.array([SW.W_func(vhp, vhm, theta, z) for z in zg]))
    med = np.median(Wz)
    zs = []
    for i in range(1, ng - 1):
        if Wz[i] < Wz[i-1] and Wz[i] < Wz[i+1] and Wz[i] < 0.02 * med:
            # parabolic refine
            y0,y1,y2 = Wz[i-1],Wz[i],Wz[i+1]
            dz = 0.5*(y0-y2)/(y0-2*y1+y2+1e-30)*(zg[1]-zg[0])
            zs.append(zg[i]+dz)
    return zs, lam_a

def track(branches, zs):
    """continuity matching: assign zeros to existing branches by nearest, else new branch."""
    if not branches:
        return [[z] for z in zs]
    last = [next((v for v in reversed(b) if v is not None), 1e9) for b in branches]
    used = set()
    for b in branches: b.append(None)
    for z in sorted(zs):
        j = min(range(len(last)), key=lambda j: abs(z-last[j]) if j not in used else 9e9)
        if abs(z-last[j]) < 1.5 and j not in used:
            branches[j][-1] = z; used.add(j)
        else:
            branches.append([None]*(len(branches[0])-1)+[z])
    return branches

if __name__ == "__main__":
    avals = sorted(set(
        list(np.arange(0.30, 1.02, 0.02)) +
        [t+d for t in THRESH for d in (-0.01,-0.005,0.0,0.005,0.01)]))
    avals = [round(a,4) for a in avals if 0.29 <= a <= 1.01]
    out = {}
    for theta in [0.0, np.pi/2, np.pi, 1.5*np.pi]:
        print(f"\n##### theta = {theta:.3f} #####", flush=True)
        branches = []
        rows = []
        for a in avals:
            try:
                zs, lam_a = zeros_at(a, theta)
            except Exception as e:
                print(f"  a={a}: FAILED {e}"); continue
            zs10 = [z for z in zs if z < 58][:14]
            branches = track(branches, zs10)
            gaps = np.diff(zs10)
            mingap = float(gaps.min()) if len(gaps) else float('nan')
            near = min((abs(a-t) for t in THRESH), default=9)
            mark = " <-- THRESHOLD q=%d" % THRESH[min(THRESH,key=lambda t:abs(a-t))] if near < 0.006 else ""
            print(f"  a={a:6.3f} lam_a={lam_a:+.2e} n={len(zs10):2d} mingap={mingap:5.2f} "
                  f"z1..4={[round(z,2) for z in zs10[:4]]}{mark}", flush=True)
            rows.append({"a": a, "lam_a": float(lam_a), "zeros": [float(z) for z in zs10],
                         "mingap": mingap})
        out[f"{theta:.4f}"] = rows
    json.dump(out, open("results/suzuki_gap_flow.json","w"), indent=1)
    print("\n-> results/suzuki_gap_flow.json")
    # summary: gap trajectories at theta=0 vs zeta gaps
    print("\nSUMMARY theta=0: last-a gaps vs zeta gaps", flush=True)
    rows = out["0.0000"]
    if rows:
        zl = rows[-1]["zeros"]; gl = list(np.round(np.diff(zl),2))
        zg = list(np.round(np.diff(GAMMAS),2))
        print(f"  final a={rows[-1]['a']}: gaps {gl[:9]}")
        print(f"  zeta gaps:                {zg[:9]}")

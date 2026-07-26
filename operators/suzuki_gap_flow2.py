#!/usr/bin/env python3
"""suzuki_gap_flow2.py -- §242: plateau-velocity scaling v_j(a) at each zeta zero, a up to 2.6.

Uses the exact-Q pipeline (suzuki_exact_Q, validated to 1e-13 vs 1-D arbiter; lambda_a = 0 to
machine floor for a>=0.8 = the critical PSD knife-edge). Observable per (theta, gamma_j):
  z_near(a) = W-zero nearest gamma_j;  v_j(a) = d z_near/da (finite difference).
Parked branch -> |v_j| small; handoff -> spike (nearest-zero identity switches). The §241
discriminator: RH predicts min-velocity in successive a-windows DECAYS (real resonance, width 0);
an off-line zero of width delta would floor at ~delta. Also records nearest-zero OFFSET
|z_near - gamma_j| per a (the §45b convergence, now to a=2.6 instead of 2.0).
"""
import numpy as np, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_exact_Q as SQ

GAMMAS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719]
ZLO, ZHI, NZ = 10.0, 44.0, 3600

def wzeros(a, theta, K, NG=320):
    W, lam_a = SQ.W_family(a, K=K, NG=NG)
    zg = np.linspace(ZLO, ZHI, NZ)
    Wz = np.abs(np.array([W(theta, z) for z in zg]))
    med = np.median(Wz)
    zs = []
    for i in range(1, NZ-1):
        if Wz[i] < Wz[i-1] and Wz[i] < Wz[i+1] and Wz[i] < 0.05*med:
            y0, y1, y2 = Wz[i-1], Wz[i], Wz[i+1]
            dz = 0.5*(y0-y2)/(y0-2*y1+y2+1e-30)*(zg[1]-zg[0])
            zs.append(zg[i]+dz)
    return zs

if __name__ == "__main__":
    avals = [round(x, 3) for x in np.arange(0.90, 2.61, 0.01)]
    out = {}
    for theta in [0.0, np.pi/2, np.pi, 1.5*np.pi]:
        rows = []
        print(f"##### theta={theta:.3f} #####", flush=True)
        for a in avals:
            K = int(np.ceil((ZHI+6)*2*a/np.pi)) + 40
            try:
                zs = wzeros(a, theta, K)
            except Exception as e:
                print(f" a={a} FAILED {e}", flush=True); continue
            near = {}
            for j, gm in enumerate(GAMMAS):
                if zs:
                    zn = min(zs, key=lambda z: abs(z-gm))
                    near[j] = zn
            rows.append({"a": a, "K": K, "n": len(zs), "near": {str(j): float(v) for j, v in near.items()}})
            if abs(a*100 - round(a*100)) < 1e-9 and (round(a*100) % 10 == 0):
                offs = [f"g{j}:{near[j]-GAMMAS[j]:+.3f}" for j in range(4) if j in near]
                print(f" a={a:5.2f} K={K:3d} nz={len(zs):2d}  {' '.join(offs)}", flush=True)
        out[f"{theta:.4f}"] = rows
    json.dump(out, open("../results/suzuki_gap_flow2.json", "w"))
    print("-> results/suzuki_gap_flow2.json", flush=True)

    # velocity summary per zero per theta: min |dz/da| over sliding windows of 0.15 in a
    print("\nVELOCITY SUMMARY (min |v| per a-window; RH => decays with a)", flush=True)
    for th, rows in out.items():
        print(f" theta={th}:")
        for j, gm in enumerate(GAMMAS[:5]):
            pts = [(r["a"], r["near"][str(j)]) for r in rows if str(j) in r["near"]]
            vs = []
            for i in range(1, len(pts)):
                da = pts[i][0]-pts[i-1][0]
                if 0 < da < 0.02:
                    vs.append((pts[i][0], abs((pts[i][1]-pts[i-1][1])/da)))
            if not vs: continue
            wins = {}
            for a, v in vs:
                w = round(np.floor(a/0.25)*0.25, 2)
                wins.setdefault(w, []).append(v)
            minv = {w: round(min(vv), 4) for w, vv in sorted(wins.items())}
            print(f"  gamma_{j+1}={gm}: min|v| per 0.25-window: {minv}")

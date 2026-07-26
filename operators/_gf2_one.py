import sys, json, numpy as np, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_gap_flow2 as GF
ti = int(sys.argv[1]); theta = [0.0, np.pi/2, np.pi, 1.5*np.pi][ti]
avals = [round(x,3) for x in np.arange(0.90, 2.61, 0.01)]
rows = []
for a in avals:
    K = int(np.ceil((GF.ZHI+6)*2*a/np.pi)) + 40
    try: zs = GF.wzeros(a, theta, K)
    except Exception as e: print(f"a={a} FAIL {e}", flush=True); continue
    near = {str(j): float(min(zs, key=lambda z: abs(z-gm))) for j, gm in enumerate(GF.GAMMAS) if zs}
    rows.append({"a": a, "K": K, "n": len(zs), "near": near})
    if round(a*100) % 20 == 0: print(f"th{ti} a={a:.2f} nz={len(zs)}", flush=True)
    json.dump(rows, open(f"../results/gf2_theta{ti}.json", "w"))
print(f"th{ti} DONE", flush=True)

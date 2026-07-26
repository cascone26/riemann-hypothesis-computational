import sys, json, numpy as np, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open("_gf2_fast.py").read().split('if __name__')[0])
lo, hi, tag = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
avals=[round(x,4) for x in np.arange(lo, hi, 0.005)]
rows=[]; t0=time.time()
for a in avals:
    K, out = scan_a(a)
    rows.append({"a": a, "K": K, "out": {str(k): v for k,v in out.items()}})
    json.dump(rows, open(f"../results/gf2ext_{tag}.json","w"))
print(f"{tag} DONE {len(rows)} in {time.time()-t0:.0f}s", flush=True)

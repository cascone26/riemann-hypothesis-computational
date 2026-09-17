import sys, time
sys.argv = ['x']  # avoid argparse
from arb_ldlt_margin import certify_margin
print("=== c-direction margin at HIGH precision (N=64) — crux of c->inf limit ===", flush=True)
print(f"{'c':>5} {'prec':>6} {'n_neg':>6} {'undet':>7} {'min_pivot_lower':>16} {'status'}", flush=True)
prev=None
for c in [14, 25, 50, 100]:
    for prec in [3000]:
        t0=time.time()
        r = certify_margin(c, 64, prec)
        ml = r["min_pivot_lower"]
        status = "PD" if r["certified_pd"] else (f"UNDET@{r['undetermined_pivot']}" if r['undetermined_pivot'] is not None else f"INDEF n_neg={r['n_neg']}")
        print(f"{c:5d} {prec:6d} {r['n_pos']-r['dim']+r['n_pos'] if False else r['n_neg']:6d} {str(r['undetermined_pivot']):>7} {('%.4e'%ml) if ml is not None else 'n/a':>16} {status}  ({time.time()-t0:.0f}s)", flush=True)

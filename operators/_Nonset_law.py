import numpy as np, math
import operators.connes_completeness_highT as C

def min_N(T, lsq, Ngrid):
    lo, hi = T-6, T+6
    for N in Ngrid:
        r = C.run(lsq, N, lo, hi)
        if r['ratio'] >= 0.999:
            return N, r['secs']
    return None, None

# fix generous prime reach; find minimum dimension N for full fill-in at each height T
LSQ = 150
Ts = [25, 40, 55, 70, 90]
print(f"prime reach LSQ={LSQ} fixed;  min N (dimension) for ratio=1 vs height T")
print(f"{'T':>4} {'minN':>5} {'N/T':>6}")
pts=[]
for T in Ts:
    # candidate N grid around the ~T*log(LSQ)/2pi expectation
    base = T*math.log(LSQ)/(2*math.pi)
    grid = sorted(set(int(base*f) for f in [0.4,0.55,0.7,0.85,1.0,1.15,1.3,1.5]))
    grid = [n for n in grid if n>=4]
    N,_ = min_N(T, LSQ, grid)
    if N: pts.append((T,N)); print(f"{T:>4} {N:>5} {N/T:>6.2f}", flush=True)
    else: print(f"{T:>4}  none  (grid {grid})", flush=True)

if len(pts)>=2:
    Tarr=np.array([p[0] for p in pts]); Narr=np.array([p[1] for p in pts])
    a,b=np.polyfit(Tarr,Narr,1)
    print(f"\nLINEAR FIT: minN = {a:.3f}*T + {b:.2f}   (dimension cost LINEAR in T)")
    lo=math.log; la,lb=np.polyfit(np.log(Tarr),np.log(Narr),1)
    print(f"LOGLOG: minN ~ T^{la:.2f}  (exponent; ~1 => linear, polynomial completeness cost)")

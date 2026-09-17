import numpy as np, math
import operators.connes_completeness_highT as C

def Nfor(T, lsq):
    L = math.log(lsq)
    return int(math.ceil(1.35 * T * L / (2*math.pi))) + 8

def onset(T, lsq_ladder):
    lo, hi = T-6, T+6
    row = []
    onset_lsq = None
    for lsq in lsq_ladder:
        N = Nfor(T, lsq)
        r = C.run(lsq, N, lo, hi)
        row.append((lsq, N, r['ratio'], r['secs']))
        if r['ratio'] >= 0.999 and onset_lsq is None:
            onset_lsq = lsq
            break
    return onset_lsq, row

Ts       = [25,   40,   55,    70,    90]
ladders  = [
    [8,14,25],
    [14,25,50,100],
    [25,50,100,200,400],
    [50,100,200,400,800],
    [100,300,700,1500],
]
print(f"{'T':>4} {'onsetLSQ':>9} {'log(onset)':>10}  detail(lsq:ratio)")
pts=[]
for T,lad in zip(Ts,ladders):
    o,row = onset(T,lad)
    detail = " ".join(f"{l}:{rt:.2f}" for (l,_,rt,_) in row)
    if o: pts.append((T, math.log(o)))
    print(f"{T:>4} {str(o):>9} {math.log(o) if o else float('nan'):>10.3f}  {detail}", flush=True)

if len(pts)>=2:
    Tarr=np.array([p[0] for p in pts]); yarr=np.array([p[1] for p in pts])
    a,b=np.polyfit(Tarr,yarr,1)
    print(f"\nFIT: log(LSQ_onset) = {a:.4f}*T + {b:.3f}  (prime-reach ~ exp({a:.4f}*T))")
    print(f"     => primes needed up to ~{a:.4f} per unit T in log; §968 4W_onset~T/7 cross-check")

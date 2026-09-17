import numpy as np, math
import operators.connes_completeness_highT as C
# generous dimension N=1.3*T; find MIN prime reach LSQ for full fill-in vs T
def min_lsq(T, Ngen, ladder):
    lo,hi=T-6,T+6
    for lsq in ladder:
        r=C.run(lsq, Ngen, lo, hi)
        if r['ratio']>=0.999: return lsq
    return None
Ts=[30,50,70,90]
print("generous N=ceil(1.3*T*log(120)/2pi); MIN prime reach LSQ for ratio=1 vs T")
print(f"{'T':>4} {'N':>4} {'minLSQ':>7}")
for T in Ts:
    Ngen=int(math.ceil(1.3*T*math.log(120)/(2*math.pi)))+6
    ladder=[6,10,14,25,50,100,200,400]
    m=min_lsq(T,Ngen,ladder)
    print(f"{T:>4} {Ngen:>4} {str(m):>7}", flush=True)

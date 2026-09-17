import numpy as np, math
import operators.connes_completeness_highT as C
from operators._snap_decay import maxerr

def decay_rate(Tc, LSQ=250, half=4):
    lo,hi=Tc-half,Tc+half
    thr=0.925*Tc
    Ns=sorted(set(int(thr*f) for f in [0.6,0.7,0.8,0.9,1.0,1.1,1.25]))
    d=[]
    for N in Ns:
        e,_=maxerr(LSQ,N,lo,hi)
        d.append((N,e))
    tr=[(N,e) for (N,e) in d if 1e-10<e<1.0]
    if len(tr)>=3:
        Na=np.array([x[0] for x in tr]); ea=np.log([x[1] for x in tr])
        a,b=np.polyfit(Na,ea,1); return a, d
    return None, d

print(f"{'T':>4} {'thrN':>5} {'decay_rate(/dim)':>16}  detail(N:err)")
for Tc in [43, 70, 95]:
    a,d=decay_rate(Tc)
    det=" ".join(f"{N}:{e:.0e}" for (N,e) in d)
    print(f"{Tc:>4} {0.925*Tc:>5.0f} {(-a if a else float('nan')):>16.3f}  {det}", flush=True)

import numpy as np, math
import operators.connes_completeness_highT as C

def maxerr(lsq, N, lo, hi):
    """max snap-error over ALL true zeros in [lo,hi] (nearest got), regardless of tol."""
    L=float(np.log(lsq)); VM=C.vm_table(lsq)
    tau,idx=C.build_tau(N,L,VM)
    ev,evec=np.linalg.eigh(tau)
    k0=int(np.argmin(np.abs(ev))); D=2*N+1; v=evec[:,k0]
    even=np.array([(v[i]+v[D-1-i])/2 for i in range(D)])
    s=even.sum()/np.sqrt(L)
    if abs(s)>1e-30: even=even/s
    xh=C.xi_hat_factory(even,idx,L)
    got=C.find_zeros(xh,lo,hi,max(4000,int((hi-lo)*120)))
    tin=C.ZTRUE[(C.ZTRUE>=lo)&(C.ZTRUE<=hi)]
    errs=[min(abs(g-tz) for g in got) if got else 9.9 for tz in tin]
    return max(errs) if errs else float('nan'), float(ev[0])

# window around T~43 (2 zeros: 40.92, 43.33); threshold N ~ 0.925*43 ~ 40
lo,hi=38,46
LSQ=200
print(f"window [{lo},{hi}], LSQ={LSQ}: max snap-error vs dimension N (threshold ~{0.925*((lo+hi)/2):.0f})")
print(f"{'N':>4} {'maxErr':>10} {'eps_min':>11}")
data=[]
for N in [16,20,24,28,32,36,40,44,48,56]:
    e,eps=maxerr(LSQ,N,lo,hi)
    data.append((N,e)); print(f"{N:>4} {e:>10.3e} {eps:>11.3e}", flush=True)

# fit exponential decay in the transition (before float floor ~1e-11)
d=[(N,e) for (N,e) in data if 1e-10 < e < 1.0]
if len(d)>=3:
    Narr=np.array([x[0] for x in d]); ee=np.log(np.array([x[1] for x in d]))
    a,b=np.polyfit(Narr,ee,1)
    print(f"\nTRANSITION FIT: log(maxErr) = {a:.3f}*N + {b:.2f}  => error ~ exp({a:.3f}*N)")
    print(f"  {'GEOMETRIC/EXPONENTIAL decay in N' if a<-0.1 else 'slower-than-exponential'} — decay rate {-a:.3f}/dim")

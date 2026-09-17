import sys, math
sys.argv=['x']
sys.path.insert(0,"../external/groskin_2607.02828/code")
from flint import ctx
from arb_ldlt_certify import build_arb_tau
from arb_ldlt_margin import inertia_with_margin
print("=== alpha(N): c-decay exponent vs N (prec=2500, c=25,50,100) ===",flush=True)
def relmarg(c,N,prec):
    ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec)
    npos,nneg,undet,mlow,mmid=inertia_with_margin(A,DIM)
    if mlow is None or nneg>0: return None
    maxd=max(abs(float(A[i,i].mid())) for i in range(DIM))
    return mlow/maxd
alphas=[]
for N in [32,48,64,80]:
    rows=[(c,relmarg(c,N,2500)) for c in [25,50,100]]
    valid=[(math.log10(c),math.log10(r)) for c,r in rows if r and r>0]
    if len(valid)>=2:
        (x0,y0),(x1,y1)=valid[0],valid[-1]; slope=(y1-y0)/(x1-x0)
    else: slope=float('nan')
    alphas.append((N,slope))
    print(f"N={N:3d}: alpha={slope:.2f}  "+"  ".join(f"c{c}:{('%.1e'%r) if r else 'na'}" for c,r in rows),flush=True)
# fit alpha vs N
va=[(N,a) for N,a in alphas if not math.isnan(a)]
if len(va)>=2:
    import numpy as np
    N_arr=np.array([x[0] for x in va]); a_arr=np.array([-x[1] for x in va])
    # linear fit |alpha| = p*N + q
    A2=np.vstack([N_arr,np.ones_like(N_arr)]).T
    (p,q),*_=np.linalg.lstsq(A2,a_arr,rcond=None)
    # power fit log|alpha| = k log N + b
    (k,b),*_=np.linalg.lstsq(np.vstack([np.log(N_arr),np.ones_like(N_arr)]).T,np.log(a_arr),rcond=None)
    print(f"\n|alpha| vs N:  linear |a|={p:.3f}*N+{q:.2f}  |  power |a|~N^{k:.2f}",flush=True)
    print(f"Interpretation: linear-in-N ⟹ joint margin ~ c^(-p*N), i.e. log-margin ~ -p*N*log c (product law).",flush=True)

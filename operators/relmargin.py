import sys, math
sys.argv=['x']
sys.path.insert(0,"../external/groskin_2607.02828/code")
from flint import ctx
from arb_ldlt_certify import build_arb_tau
from arb_ldlt_margin import inertia_with_margin
print("=== RELATIVE positivity margin vs c (N=64, prec=3000): conditioning or danger? ===",flush=True)
print(f"{'c':>5} {'min_pivot':>13} {'max|diag|':>13} {'min/maxdiag':>13} {'log10(rel)':>11}",flush=True)
for c in [14,25,50,100]:
    ctx.prec=3000
    A,DIM = build_arb_tau(c,64,3000)
    npos,nneg,undet,mlow,mmid = inertia_with_margin(A,DIM)
    maxd=0.0
    for i in range(DIM):
        v=abs(float(A[i,i].mid()))
        if v>maxd: maxd=v
    rel = mlow/maxd if (mlow is not None and maxd>0) else None
    lr = math.log10(rel) if (rel and rel>0) else float('nan')
    print(f"{c:5d} {mlow:13.3e} {maxd:13.3e} {rel:13.3e} {lr:11.3f}",flush=True)

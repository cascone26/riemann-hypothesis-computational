import sys, math
sys.argv=['x']
sys.path.insert(0,"../external/groskin_2607.02828/code")
from flint import ctx
from arb_ldlt_certify import build_arb_tau
from arb_ldlt_margin import inertia_with_margin
print("=== EXTENDED c-scan (N=64, prec=6000): does deceleration continue, margin stay positive? ===",flush=True)
print(f"{'c':>5} {'min_pivot':>13} {'max|diag|':>13} {'min/maxdiag':>13} {'log10(rel)':>11} {'PD?'}",flush=True)
for c in [100,150,200,300]:
    ctx.prec=6000
    A,DIM = build_arb_tau(c,64,6000)
    npos,nneg,undet,mlow,mmid = inertia_with_margin(A,DIM)
    maxd=max(abs(float(A[i,i].mid())) for i in range(DIM))
    rel = mlow/maxd if (mlow is not None and maxd>0) else None
    lr = math.log10(rel) if (rel and rel>0) else float('nan')
    pd = "PD" if (nneg==0 and undet is None) else (f"UNDET@{undet}" if undet is not None else f"INDEF{nneg}")
    print(f"{c:5d} {(mlow if mlow else 0):13.3e} {maxd:13.3e} {(rel if rel else 0):13.3e} {lr:11.3f} {pd}",flush=True)

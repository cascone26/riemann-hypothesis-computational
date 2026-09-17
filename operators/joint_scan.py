import sys, math
sys.argv=['x']
sys.path.insert(0,"../external/groskin_2607.02828/code")
from flint import ctx
from arb_ldlt_certify import build_arb_tau
from arb_ldlt_margin import inertia_with_margin
print("=== JOINT (N,c) margin: does c-decay exponent depend on N? ===",flush=True)
def relmarg(c,N,prec):
    ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec)
    npos,nneg,undet,mlow,mmid=inertia_with_margin(A,DIM)
    if mlow is None or nneg>0: return None
    maxd=max(abs(float(A[i,i].mid())) for i in range(DIM))
    return mlow/maxd
for N in [32,64,96]:
    rows=[]
    for c in [25,50,100]:
        r=relmarg(c,N,2000)
        rows.append((c,r))
    # fit log-log slope over c=25->100
    valid=[(math.log10(c),math.log10(r)) for c,r in rows if r and r>0]
    if len(valid)>=2:
        (x0,y0),(x1,y1)=valid[0],valid[-1]
        slope=(y1-y0)/(x1-x0)
    else: slope=float('nan')
    print(f"N={N:3d}: " + "  ".join(f"c={c}:{('%.2e'%r) if r else 'n/a'}" for c,r in rows) + f"   | c-decay exponent alpha={slope:.2f}",flush=True)
print("\nIf alpha grows with N: joint limit subtle (margin decays faster at higher band). If ~constant: cleaner.",flush=True)

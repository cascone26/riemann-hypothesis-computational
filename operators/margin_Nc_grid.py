#!/usr/bin/env python3
"""Joint (N,c) law of the certified positivity margin of tau.

Established: at fixed N=48 the certified min-pivot margin obeys margin(c) ~ c^(-p) with p~32 (log10-log10,
R^2=0.998). Question: is the exponent p N-dependent? Map margin over an (N,c) grid, fit p(N). Together with
the earlier N-plateau result this pins how positivity behaves in the joint N,c -> infinity limit that IS RH.
"""
import sys, os, time, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, ctx
from arb_ldlt_certify import build_arb_tau

def margin(c, N, prec):
    ctx.prec = prec
    A, DIM = build_arb_tau(c, N, prec)
    d=[None]*DIM; Lf=[[arb(0)]*DIM for _ in range(DIM)]; mlow=None; undet=None; nneg=0
    for i in range(DIM):
        s=A[i,i]
        for k in range(i): s=s-Lf[i][k]*Lf[i][k]*d[k]
        d[i]=s
        if s>0:
            low=float(s.mid())-float(s.rad())
            if mlow is None or low<mlow: mlow=low
        elif s<0: nneg+=1
        else: undet=i; break
        for j in range(i+1,DIM):
            t=A[j,i]
            for k in range(i): t=t-Lf[j][k]*Lf[i][k]*d[k]
            Lf[j][i]=t/d[i]
    return mlow, undet, nneg

def fit_loglog(cs, ys):
    xs=[math.log10(c) for c in cs]; ly=[math.log10(y) for y in ys]
    n=len(xs);sx=sum(xs);sy=sum(ly);sxx=sum(x*x for x in xs);sxy=sum(x*y for x,y in zip(xs,ly))
    m=(n*sxy-sx*sy)/(n*sxx-sx*sx); b=(sy-m*sx)/n
    ss=sum((y-(m*x+b))**2 for x,y in zip(xs,ly)); tot=sum((y-sy/n)**2 for y in ly)
    return m,(1-ss/tot if tot else 0)

def main():
    PREC=6000
    Ns=[24,48,96]; cs=[8,25,100,400]
    print(f"=== (N,c) margin grid, prec={PREC} bits ===",flush=True)
    grid={}
    for N in Ns:
        ys=[]; row_cs=[]
        for c in cs:
            t0=time.time(); mlow,undet,nneg=margin(c,N,PREC)
            ok = undet is None and nneg==0 and mlow and mlow>0
            print(f"  N={N:3d} c={c:4d}  margin={('%.3e'%mlow) if mlow else 'n/a':>11}  "
                  f"log10={('%.2f'%math.log10(mlow)) if ok else 'n/a':>7}  {'PD' if ok else 'UNDET/INDEF'}  ({time.time()-t0:.0f}s)",flush=True)
            if ok: ys.append(mlow); row_cs.append(c)
        if len(ys)>=3:
            p,r2=fit_loglog(row_cs,ys)
            grid[N]=dict(exponent=p,R2=r2,cs=row_cs,margins=ys)
            print(f"    -> N={N}: margin ~ c^({p:+.2f})  R^2={r2:.3f}",flush=True)
    json.dump(grid,open("results/margin_Nc_grid.json","w"),indent=2)
    print("\n=== EXPONENT vs N ===",flush=True)
    for N in sorted(grid): print(f"  N={N:3d}  exponent p={grid[N]['exponent']:+.2f}  R^2={grid[N]['R2']:.3f}",flush=True)
    if len(grid)>=2:
        Nv=sorted(grid); ps=[grid[N]['exponent'] for N in Nv]
        spread=max(ps)-min(ps)
        print(f"\n  exponent spread across N: {spread:.2f}  =>  "
              f"{'N-INDEPENDENT (a universal prime-side collapse law)' if abs(spread)<3 else 'N-DEPENDENT (joint law couples N and c)'}",flush=True)
    print("wrote results/margin_Nc_grid.json",flush=True)

if __name__=="__main__":
    main()

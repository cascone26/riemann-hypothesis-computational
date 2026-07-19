#!/usr/bin/env python3
"""How fast does the certified positivity margin of tau collapse as the prime cutoff c grows?

Finding from arb_ldlt_margin.py: at fixed N the min LDL^T pivot (a rigorous positivity margin) drops
sharply with c and 1000 bits stops certifying at c>=25. Here we use HIGH, fixed precision to actually
resolve the margin for a c-ladder at fixed N, and fit the decay law margin(c). This localizes the real
crux (the prime side) quantitatively — a genuine measurement, not a re-verification.
"""
import sys, os, time, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, ctx
from arb_ldlt_certify import build_arb_tau

def margin(c, N, prec):
    ctx.prec = prec
    A, DIM = build_arb_tau(c, N, prec)
    d = [None]*DIM; Lf=[[arb(0)]*DIM for _ in range(DIM)]
    mlow=None; undet=None; nneg=0
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
    return dict(c=c,N=N,dim=DIM,prec=prec,min_pivot_lower=mlow,undetermined_pivot=undet,n_neg=nneg,
                certified_pd=bool(undet is None and nneg==0))

def main():
    N=48; PREC=6000    # dim=97; 6000 bits (~1800 digits) — enough to resolve deep margins
    cs=[8,14,25,50,100,200,400]
    print(f"=== certified margin(c) at N={N}, prec={PREC} bits ===",flush=True)
    rows=[]
    for c in cs:
        t0=time.time(); r=margin(c,N,PREC); rows.append(r)
        tag="PD" if r["certified_pd"] else (f"UNDET@{r['undetermined_pivot']}" if r["undetermined_pivot"] is not None else f"INDEF n_neg={r['n_neg']}")
        mp=r["min_pivot_lower"]
        print(f"  c={c:4d}  min_pivot_lower={('%.4e'%mp) if mp else 'n/a':>12}  log10={('%.2f'%math.log10(mp)) if mp and mp>0 else 'n/a':>7}  [{tag}]  ({time.time()-t0:.0f}s)",flush=True)
    json.dump(rows,open("results/margin_c_decay.json","w"),indent=2)
    # fit log10(margin) vs c  and vs log(c)  and vs pi(c)-ish (number of primes)
    good=[(r["c"],r["min_pivot_lower"]) for r in rows if r["certified_pd"] and r["min_pivot_lower"] and r["min_pivot_lower"]>0]
    if len(good)>=3:
        cs2=[c for c,_ in good]; ly=[math.log10(y) for _,y in good]
        def lin(xs):
            n=len(xs);sx=sum(xs);sy=sum(ly);sxx=sum(x*x for x in xs);sxy=sum(x*y for x,y in zip(xs,ly))
            m=(n*sxy-sx*sy)/(n*sxx-sx*sx); b=(sy-m*sx)/n
            ss=sum((y-(m*x+b))**2 for x,y in zip(xs,ly)); tot=sum((y-sy/n)**2 for y in ly)
            return m,b,(1-ss/tot if tot else 0)
        mc,_,r2c=lin(cs2)
        mlc,_,r2lc=lin([math.log(c) for c in cs2])
        print("\n=== DECAY LAW (log10 margin vs ...) ===",flush=True)
        print(f"  vs c        : slope={mc:+.4f}/unit-c   R^2={r2c:.3f}   (=> margin ~ 10^({mc:+.3f} c), exponential in c)",flush=True)
        print(f"  vs log c    : slope={mlc:+.3f}/e-fold  R^2={r2lc:.3f}   (=> margin ~ c^({mlc*math.log(10)/math.log(10):+.2f}) power-law)",flush=True)
        print("  Better R^2 names the true collapse law of the prime-side positivity margin.",flush=True)
    print("wrote results/margin_c_decay.json",flush=True)

if __name__=="__main__":
    main()

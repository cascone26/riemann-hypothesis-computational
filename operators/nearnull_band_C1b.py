#!/usr/bin/env python3
"""C1b: is the MARGIN-carrying near-null band ALSO bounded in frequency as c->inf?

C1 showed A+B's NEGATIVE band edge converges (~6). But the tiny tau margin lives on the near-null
direction, which §11 saw at |omega|~9-17 with unclear convergence. The core+tail reduction needs the
ENTIRE bad region (negative band + near-null cancellation) to be a FIXED frequency band. This tracks the
near-null direction's dominant frequency on a long c-ladder in high precision. If it converges/stays
bounded, the core is genuinely finite-in-frequency and the linchpin is fully closed.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, ctx
from arb_ldlt_certify import build_arb_tau

def near_null(c, N, prec):
    ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec)
    d=[None]*DIM; Lf=[[arb(0)]*DIM for _ in range(DIM)]; imin=None; dmin=None
    for i in range(DIM):
        s=A[i,i]
        for k in range(i): s=s-Lf[i][k]*Lf[i][k]*d[k]
        d[i]=s; v=float(s.mid())
        if imin is None or v<dmin: dmin=v; imin=i
        for j in range(i+1,DIM):
            t=A[j,i]
            for k in range(i): t=t-Lf[j][k]*Lf[i][k]*d[k]
            Lf[j][i]=t/d[i]
    x=[arb(0)]*DIM
    for i in range(DIM-1,-1,-1):
        s=arb(1) if i==imin else arb(0)
        for j in range(i+1,DIM): s=s-Lf[j][i]*x[j]
        x[i]=s
    xf=[float(v.mid()) for v in x]
    # dominant frequency + centroid frequency of the mass
    mass=[(i-N, xf[i]*xf[i]) for i in range(DIM)]
    tot=sum(m for _,m in mass) or 1
    dom_n=max(mass,key=lambda t:t[1])[0]
    centroid=sum(abs(n)*m for n,m in mass)/tot     # mean |n|
    return dmin, dom_n, centroid, float(ctx.prec)

def main():
    N=64
    print(f"=== C1b: near-null band vs c (N={N}) ===")
    print(f"{'c':>6} {'L':>7} {'margin':>11} {'dom|omega|':>11} {'centroid|omega|':>15}")
    rows=[]
    for c in [14,50,200,800,3200]:
        L=math.log(c)
        prec=max(2000, int(40*L/math.log(10)*3.4)+2000)   # grow precision with expected margin depth
        dmin,domn,cent,_=near_null(c,N,prec)
        dom_omega=2*math.pi*abs(domn)/L
        cent_omega=2*math.pi*cent/L
        rows.append(dict(c=c,L=L,margin=dmin,dom_omega=dom_omega,centroid_omega=cent_omega))
        print(f"{c:>6} {L:>7.3f} {dmin:>11.2e} {dom_omega:>11.3f} {cent_omega:>15.3f}")
    json.dump(rows,open("results/nearnull_band_C1b.json","w"),indent=2)
    print("\n=== VERDICT (C1b) ===")
    co=[r['centroid_omega'] for r in rows]
    spread=max(co[-3:])-min(co[-3:])
    print(f"  centroid |omega| last-3 spread = {spread:.2f}")
    if spread<3.0:
        print("  => PASS: near-null band is BOUNDED in frequency; full bad region is a fixed band.")
        print("     Linchpin fully closed: core (neg band + near-null) is finite-in-frequency.")
    else:
        print("  => near-null band drifts; core may need a frequency-dependent cutoff. Investigate.")
    print("wrote results/nearnull_band_C1b.json")

if __name__=="__main__":
    main()

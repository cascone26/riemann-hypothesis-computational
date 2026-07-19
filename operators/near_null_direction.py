#!/usr/bin/env python3
"""WHAT nearly breaks positivity? Extract the near-null direction of the real Weil operator tau.

From tau = L D L^T (Arb interval LDL^T), the vector achieving the SMALLEST quadratic form is
x = L^{-T} e_{i*}, where i* is the index of the smallest pivot d[i*]; then x^T tau x = d[i*] (tiny).
This x is the direction that most nearly makes the Weil positivity form vanish — the numerical stand-in
for "the mode the operator wants as a null vector." As the prime cutoff c grows the margin collapses
(LEDGER §10), so this direction is where the action is. We extract it exactly, then characterize it:

  - which basis frequencies n (freq omega = 2*pi*n / log c) carry its mass?
  - is it localized (a few modes) or spread?
  - do its dominant frequencies track anything recognizable across c?

Honest: this is the min-PIVOT direction (exact, from the certified factorization), a faithful proxy for
the near-null eigenvector. Prints the mass profile + dominant modes for a c-ladder. No claim beyond what
the numbers show.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, ctx
from arb_ldlt_certify import build_arb_tau

def ldlt(A, DIM):
    d=[None]*DIM; Lf=[[arb(0)]*DIM for _ in range(DIM)]
    imin=None; dmin=None
    for i in range(DIM):
        s=A[i,i]
        for k in range(i): s=s-Lf[i][k]*Lf[i][k]*d[k]
        d[i]=s
        val=float(s.mid())
        if imin is None or val<dmin: dmin=val; imin=i
        for j in range(i+1,DIM):
            t=A[j,i]
            for k in range(i): t=t-Lf[j][k]*Lf[i][k]*d[k]
            Lf[j][i]=t/d[i]
    return d, Lf, imin, dmin

def near_null(c, N, prec):
    ctx.prec=prec
    A,DIM=build_arb_tau(c,N,prec)
    d,Lf,imin,dmin=ldlt(A,DIM)
    # solve L^T x = e_{imin}. L is unit-lower-triangular (Lf[j][k], j>k; diag 1).
    # L^T is unit-upper; back-substitute from bottom.
    x=[arb(0)]*DIM
    for i in range(DIM-1,-1,-1):
        rhs=arb(1) if i==imin else arb(0)
        s=rhs
        for j in range(i+1,DIM):
            s=s-Lf[j][i]*x[j]   # (L^T)_{i j} = L_{j i} = Lf[j][i]
        x[i]=s                  # diag of L^T is 1
    xf=[float(v.mid()) for v in x]
    norm=math.sqrt(sum(v*v for v in xf)) or 1.0
    xf=[v/norm for v in xf]
    return xf, imin, dmin, DIM

def analyze(c, N=64, prec=4000, zeros=(14.1347,21.0220,25.0109,30.4249,32.9351)):
    xf, imin, dmin, DIM = near_null(c, N, prec)
    L=math.log(c)
    # basis index i -> n = i-N -> frequency omega = 2*pi*n / L
    mass=[(i-N, xf[i]*xf[i]) for i in range(DIM)]
    mass.sort(key=lambda t:-t[1])
    top=mass[:8]
    total=sum(m for _,m in mass)
    top5frac=sum(m for _,m in mass[:5])/total if total else 0
    print(f"\nc={c:4d} (L=log c={L:.3f})  min-pivot idx i*={imin}/{DIM} (n*={imin-N})  margin~{dmin:.2e}")
    print(f"  concentration: top-5 modes hold {100*top5frac:.1f}% of the mass")
    print(f"  dominant modes (n, freq omega=2pi n/L, mass%):")
    for n,m in top:
        omega=2*math.pi*n/L
        # nearest zeta zero to |omega|
        znear=min(zeros, key=lambda g:abs(g-abs(omega)))
        tag=f"  ~zero {znear}?" if abs(abs(omega)-znear)<0.5 else ""
        print(f"    n={n:+4d}  omega={omega:+8.3f}  mass={100*m/total:5.1f}%{tag}")
    return dict(c=c, N=N, imin=imin, nstar=imin-N, margin=dmin, top5frac=top5frac,
                top=[(n,m/total) for n,m in top])

def main():
    print("=== near-null direction of tau across the prime cutoff ===")
    out=[]
    for c in [14, 25, 50, 100, 200]:
        out.append(analyze(c, N=64, prec=4000))
    # does n* (peak location) move with c in a recognizable way?
    print("\n=== n* (near-null peak index) vs c ===")
    for r in out:
        print(f"  c={r['c']:4d}  n*={r['nstar']:+4d}  top5={100*r['top5frac']:.0f}%  margin={r['margin']:.1e}")
    json.dump(out, open("results/near_null_direction.json","w"), indent=2)
    print("wrote results/near_null_direction.json")

if __name__=="__main__":
    main()

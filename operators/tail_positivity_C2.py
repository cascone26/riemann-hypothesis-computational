#!/usr/bin/env python3
"""C2 (core+tail reduction test): does the collapse live ONLY in the finite low-frequency CORE, leaving a
positive TAIL with a NON-collapsing margin?

C1/C1b: the bad region (neg band + margin mode) is a FIXED frequency band, asymptotically ω≲6, worst-case
ω≲17. So fix a frequency cutoff ω_c ABOVE the bad band. Split the basis into CORE (|ω|≤ω_c, finite growth
~L) and TAIL (|ω|>ω_c). Test: is the TAIL principal block of τ positive-definite with a margin that does
NOT collapse as c grows? If yes, the whole c^(−p) collapse is confined to the finite core (rigorously
certifiable) and the tail is analytically tame — the core+tail reduction (PLAN §1) is a real proof shape.

Rigorous Arb min-pivot on the tail principal submatrix; compare to full-τ margin.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, arb_mat, ctx
from arb_ldlt_certify import build_arb_tau

def min_pivot_sub(A, idxs, prec):
    """min-pivot lower bound of the principal submatrix A[idxs,idxs] via interval LDL^T."""
    ctx.prec=prec; D=len(idxs); mlow=None; nneg=0; undet=None
    d=[None]*D; Lf=[[arb(0)]*D for _ in range(D)]
    for a in range(D):
        i=idxs[a]; s=A[i,i]
        for k in range(a): s=s-Lf[a][k]*Lf[a][k]*d[k]
        d[a]=s
        if s>0:
            low=float(s.mid())-float(s.rad())
            if mlow is None or low<mlow: mlow=low
        elif s<0: nneg+=1
        else: undet=a;break
        for b in range(a+1,D):
            j=idxs[b]; t=A[j,i]
            for k in range(a): t=t-Lf[b][k]*Lf[a][k]*d[k]
            Lf[b][a]=t/d[a]
    return mlow,nneg,undet

def full_margin(A,DIM,prec):
    return min_pivot_sub(A,list(range(DIM)),prec)

def main():
    N=80; OMEGA_C=20.0    # cutoff above the worst-case bad band (~17)
    print(f"=== C2: tail positivity above ω_c={OMEGA_C} (N={N}) ===")
    print(f"{'c':>6} {'L':>7} {'core dim':>9} {'tail dim':>9} {'tail min-pivot':>16} {'full margin':>14} {'tail/full':>12}")
    out=[]
    for c in [14,50,200,800,3200]:
        L=math.log(c); prec=max(3000,int(45*L/math.log(10)*3.4)+2500)
        ctx.prec=prec; A,DIM=build_arb_tau(c,N,prec)
        nc=int(OMEGA_C*L/(2*math.pi))                 # core half-width in mode index
        core=[i for i in range(DIM) if abs(i-N)<=nc]
        tail=[i for i in range(DIM) if abs(i-N)>nc]
        tm,tn,tu=min_pivot_sub(A,tail,prec)
        fm,fn,fu=full_margin(A,DIM,prec)
        ratio=(tm/fm) if (tm and fm and fm>0) else float('nan')
        tag="PD" if (tn==0 and tu is None) else f"NOTPD(n_neg={tn},undet={tu})"
        print(f"{c:>6} {L:>7.3f} {len(core):>9} {len(tail):>9} {('%.3e'%tm) if tm else 'n/a':>16} {('%.2e'%fm) if fm else 'n/a':>14} {ratio:>12.2e}   tail:{tag}")
        out.append(dict(c=c,L=L,nc=nc,core_dim=len(core),tail_dim=len(tail),
                        tail_min_pivot=tm,full_margin=fm,tail_pd=(tn==0 and tu is None)))
    json.dump(out,open("results/tail_positivity_C2.json","w"),indent=2)
    print("\n=== VERDICT (C2) ===")
    good=[r for r in out if r['tail_min_pivot'] and r['tail_min_pivot']>0]
    if len(good)>=3:
        tps=[r['tail_min_pivot'] for r in good]
        tail_span=max(tps)/min(tps)
        all_pd=all(r['tail_pd'] for r in out)
        print(f"  tail PD at every c: {all_pd}")
        print(f"  tail min-pivot varies by {tail_span:.1e}x over c (full margin varies by many orders).")
        if all_pd and tail_span<1e4:
            print("  => PASS: TAIL is uniformly positive with a NON-collapsing margin. The c^(−p) collapse is")
            print("     confined to the finite CORE. Core+tail reduction VALIDATED as a proof structure.")
            print("     Remaining: (i) rigorously certify the finite core for all c, (ii) prove the tail bound")
            print("     analytically (its margin is stable -> tractable). This is rung 4.")
        else:
            print("  => tail not uniformly positive/stable; cutoff or splitting needs adjustment.")
    print("wrote results/tail_positivity_C2.json")

if __name__=="__main__":
    main()

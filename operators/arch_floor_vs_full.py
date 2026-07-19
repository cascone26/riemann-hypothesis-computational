#!/usr/bin/env python3
"""Is the archimedean side a ROBUST floor that the prime side merely grows to nearly meet?

tau = (A+B) - P.  Weil positivity <=> (A+B) >= P as forms <=> top generalized eigenvalue
lambda_max(P ; A+B) <= 1, and margin(full tau) is set by 1 - lambda_max.

Non-circular test: certify (A+B) ALONE (its own interval LDL^T) and compare its smallest pivot
(the archimedean floor mu_arch) against the smallest pivot of the FULL tau across the prime cutoff c.
  - if mu_arch stays O(1) healthy while margin(tau) collapses  => the collapse is a NEAR-CANCELLATION:
    the floor is fine, P grows to nearly meet it. The proof problem is then a clean two-sided bound:
    floor(A+B) >= something, ||P||_{A+B} <= 1 - something, uniformly in c.
  - if mu_arch ALSO collapses => the archimedean side itself is degenerating; different problem.
Rigorous Arb; prints mu_arch(c), margin_full(c), and their ratio.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, arb_mat, ctx
from arb_ldlt_certify import arb_closed_forms, arb_kappa, arb_J, prime_powers_up_to

def build_split(c, N, prec):
    ctx.prec = prec
    S, CC, XC, L = arb_closed_forms(N, c, prec)
    PI = arb.pi(); sp2 = 16*PI*PI; l2 = L*L
    pref02 = 32*L*(L/4).sinh()**2; kappa = arb_kappa(L); J = arb_J(L)
    pdata = prime_powers_up_to(c)
    weights = [arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    positions = [arb(q).log() for (q,p) in pdata]; M=len(weights)
    def Ss(nn): return S[nn] if nn>=0 else -S[-nn]
    DIM=2*N+1; Aa=arb_mat(DIM,DIM); Af=arb_mat(DIM,DIM)
    for i in range(DIM):
        n=i-N
        for j in range(i,DIM):
            m=j-N
            num=l2-sp2*m*n; den=(l2+sp2*m*m)*(l2+sp2*n*n); W02=pref02*num/den
            if n==m: WR=kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]
            else:    WR=(Ss(m)-Ss(n))/(PI*(n-m))
            Wp=arb(0)
            for idx in range(M):
                y=positions[idx]
                if n==m: q=2*(1-y/L)*(2*PI*n*y/L).cos()
                else:    q=((2*PI*m*y/L).sin()-(2*PI*n*y/L).sin())/(PI*(n-m))
                Wp+=weights[idx]*q
            arch=W02-WR; full=arch-Wp
            Aa[i,j]=arch; Aa[j,i]=arch; Af[i,j]=full; Af[j,i]=full
    return Aa, Af, DIM

def min_pivot(A, DIM):
    d=[None]*DIM; Lf=[[arb(0)]*DIM for _ in range(DIM)]; mlow=None; nneg=0; undet=None
    for i in range(DIM):
        s=A[i,i]
        for k in range(i): s=s-Lf[i][k]*Lf[i][k]*d[k]
        d[i]=s
        if s>0:
            low=float(s.mid())-float(s.rad())
            if mlow is None or low<mlow: mlow=low
        elif s<0: nneg+=1
        else: undet=i;break
        for j in range(i+1,DIM):
            t=A[j,i]
            for k in range(i): t=t-Lf[j][k]*Lf[i][k]*d[k]
            Lf[j][i]=t/d[i]
    return mlow, nneg, undet

def main():
    N=64; PREC=4000
    print(f"=== archimedean floor vs full margin (N={N}, prec={PREC}) ===")
    print(f"{'c':>5} {'mu_arch (A+B floor)':>22} {'margin_full (tau)':>20} {'floor/margin':>14}")
    out=[]
    for c in [8,14,25,50,100,200]:
        Aa,Af,DIM=build_split(c,N,PREC)
        ma,na,ua=min_pivot(Aa,DIM)
        mf,nf,uf=min_pivot(Af,DIM)
        ratio = (ma/mf) if (ma and mf and mf>0) else float('nan')
        aok="PD" if (na==0 and ua is None) else f"NOTPD(n_neg={na},undet={ua})"
        print(f"{c:>5} {('%.4e'%ma) if ma else 'n/a':>22} {('%.4e'%mf) if mf else 'n/a':>20} {ratio:>14.3e}   archA+B:{aok}")
        out.append(dict(c=c,mu_arch=ma,margin_full=mf,floor_over_margin=ratio,arch_pd=(na==0 and ua is None)))
    json.dump(out,open("results/arch_floor_vs_full.json","w"),indent=2)
    print("\n=== READ ===")
    good=[r for r in out if r['mu_arch'] and r['margin_full'] and r['mu_arch']>0 and r['margin_full']>0]
    if len(good)>=3:
        # does mu_arch stay ~flat while margin_full collapses?
        mus=[r['mu_arch'] for r in good]; mfs=[r['margin_full'] for r in good]
        arch_span=max(mus)/min(mus); full_span=max(mfs)/min(mfs)
        print(f"  archimedean floor mu_arch varies by {arch_span:.1e}x over c; full margin varies by {full_span:.1e}x.")
        if arch_span < 1e3 and full_span > 1e6:
            print("  => FLOOR IS ROBUST, full margin collapses: the positivity loss is a NEAR-CANCELLATION")
            print("     (P grows to nearly meet a healthy archimedean floor). Proof target: bound P below the")
            print("     floor uniformly in c. This is the clean, non-circular shape of the problem.")
        else:
            print("  => floor and margin move together; the archimedean side is not an independent robust floor.")
    print("wrote results/arch_floor_vs_full.json")

if __name__=="__main__":
    main()

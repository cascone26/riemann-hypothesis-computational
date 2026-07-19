#!/usr/bin/env python3
"""THE RACE: on the near-null direction, is the tiny margin a near-cancellation of a large archimedean
floor against a large prime erosion?  tau = (A+B) - P  with  (A+B)=W02-WR (archimedean+boundary),  P=Wp.

For the near-null direction x (extracted from the LDL^T smallest-pivot), compute separately
   Farch = x^T (A+B) x   (the positive floor the archimedean side provides)
   Fprime= x^T  P     x   (how much the prime sum eats)
   margin= Farch - Fprime = x^T tau x  (the tiny certified pivot-scale value)
If Farch and Fprime are both O(1) and nearly equal, positivity survives only by a delicate cancellation
whose SIZE (Farch) and CLOSENESS (margin/Farch) are the real quantities a proof must control. This is the
shape of a Weil-positivity proof: bound the floor below, bound the erosion above, show floor > erosion
uniformly in c. Here we MEASURE that race on the exact strain mode.

Shares one L=log c between the two pieces (the archimedean kernel depends on L too, so you cannot drop
primes by setting c=1). Rigorous Arb entries; reports midpoints (the ratio, not the last digit, is the point).
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, arb_mat, ctx
from arb_ldlt_certify import arb_closed_forms, arb_kappa, arb_J, prime_powers_up_to

def build_split(c, N, prec):
    """Return (A_arch, A_prime, DIM): A_arch = W02 - WR, A_prime = Wp, same L=log c."""
    ctx.prec = prec
    S, CC, XC, L = arb_closed_forms(N, c, prec)
    PI = arb.pi(); sp2 = 16*PI*PI; l2 = L*L
    pref02 = 32*L*(L/4).sinh()**2
    kappa = arb_kappa(L); J = arb_J(L)
    pdata = prime_powers_up_to(c)
    weights = [arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    positions = [arb(q).log() for (q,p) in pdata]
    M = len(weights)
    def S_signed(nn): return S[nn] if nn>=0 else -S[-nn]
    DIM = 2*N+1
    Aa = arb_mat(DIM,DIM); Ap = arb_mat(DIM,DIM)
    for i in range(DIM):
        n = i-N
        for j in range(i,DIM):
            m = j-N
            num = l2 - sp2*m*n; den = (l2+sp2*m*m)*(l2+sp2*n*n)
            W02 = pref02*num/den
            if n==m: WR = kappa + 2*CC[abs(n)] + J - (2/L)*XC[abs(n)]
            else:    WR = (S_signed(m)-S_signed(n))/(PI*(n-m))
            Wp = arb(0)
            for idx in range(M):
                y = positions[idx]
                if n==m: q = 2*(1-y/L)*(2*PI*n*y/L).cos()
                else:    q = ((2*PI*m*y/L).sin()-(2*PI*n*y/L).sin())/(PI*(n-m))
                Wp += weights[idx]*q
            arch = W02 - WR
            Aa[i,j]=arch; Aa[j,i]=arch
            Ap[i,j]=Wp;   Ap[j,i]=Wp
    return Aa, Ap, DIM

def ldlt_min_dir(A, DIM):
    """LDL^T of A(full = arch - prime); return near-null direction x (L^-T e_{i*}) as floats."""
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
    return x, imin, dmin

def quad(A, x, DIM):
    """x^T A x as an arb (x list of arb)."""
    tot=arb(0)
    for i in range(DIM):
        if float(x[i].mid())==0: continue
        row=arb(0)
        for j in range(DIM):
            row += A[i,j]*x[j]
        tot += x[i]*row
    return tot

def main():
    N=64; PREC=4000
    print(f"=== archimedean-floor vs prime-erosion RACE on the near-null mode (N={N}, prec={PREC}) ===")
    out=[]
    for c in [14,25,50,100]:
        Aa,Ap,DIM=build_split(c,N,PREC)
        # full tau = Aa - Ap ; get its near-null direction
        Afull=arb_mat(DIM,DIM)
        for i in range(DIM):
            for j in range(DIM):
                Afull[i,j]=Aa[i,j]-Ap[i,j]
        x,imin,dmin=ldlt_min_dir(Afull,DIM)
        # normalize x
        nrm=arb(0)
        for i in range(DIM): nrm+=x[i]*x[i]
        nrm=nrm.sqrt()
        x=[xi/nrm for xi in x]
        Farch=float(quad(Aa,x,DIM).mid())
        Fprime=float(quad(Ap,x,DIM).mid())
        margin=Farch-Fprime
        ratio = margin/Farch if Farch else float('nan')
        print(f"\nc={c:4d}  i*={imin}(n*={imin-N})")
        print(f"  Farch  (floor)   = {Farch:+.6e}")
        print(f"  Fprime (erosion) = {Fprime:+.6e}")
        print(f"  margin = Farch-Fprime = {margin:+.3e}   (cancellation: margin/Farch = {ratio:.2e})")
        out.append(dict(c=c,nstar=imin-N,Farch=Farch,Fprime=Fprime,margin=margin,cancel_ratio=ratio))
    json.dump(out,open("results/race_decomposition.json","w"),indent=2)
    print("\n=== READ ===")
    print("  cancel_ratio = margin/Farch  is how close the race is: 1e-3 means floor and erosion agree to")
    print("  0.1%, so positivity survives on a knife-edge whose control IS the proof problem.")
    for r in out:
        print(f"  c={r['c']:4d}  floor={r['Farch']:+.3e}  erosion={r['Fprime']:+.3e}  cancel={r['cancel_ratio']:.1e}")
    print("wrote results/race_decomposition.json")

if __name__=="__main__":
    main()

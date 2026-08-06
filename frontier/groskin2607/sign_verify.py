"""VERIFY §741: is arch+pole-only (W02-WR) positive-definite or INDEFINITE? Compute the
ALGEBRAIC min eigenvalue (with sign), not the closest-to-zero one. Compare full tau vs arch-only."""
import sys, math, numpy as np
sys.path.insert(0,'.')
from arb_ldlt_certify import arb_closed_forms, arb_J, arb_kappa, prime_powers_up_to
from flint import arb, ctx
import mpmath as mp
def build(c,N,prec,include_prime):
    ctx.prec=prec; S,CC,XC,L=arb_closed_forms(N,c,prec); PI=arb.pi(); sp2=16*PI*PI; l2=L*L
    pref02=32*L*(L/4).sinh()**2; kappa=arb_kappa(L); J=arb_J(L)
    pdata=prime_powers_up_to(c); weights=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    positions=[arb(q).log() for (q,p) in pdata]; M=len(weights)
    def Ss(nn): return S[nn] if nn>=0 else -S[-nn]
    DIM=2*N+1; mp.mp.dps=int(prec*0.29); Amp=mp.matrix(DIM,DIM)
    for i in range(DIM):
        n=i-N
        for j in range(DIM):
            m=j-N; num=l2-sp2*m*n; den=(l2+sp2*m*m)*(l2+sp2*n*n); W02=pref02*num/den
            if n==m: WR=kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]
            else: WR=(Ss(m)-Ss(n))/(PI*(n-m))
            Wp=arb(0)
            if include_prime:
                for idx in range(M):
                    y=positions[idx]
                    if n==m: qq=2*(1-y/L)*(2*PI*n*y/L).cos()
                    else: qq=((2*PI*m*y/L).sin()-(2*PI*n*y/L).sin())/(PI*(n-m))
                    Wp+=weights[idx]*qq
            Amp[i,j]=mp.mpf((W02-WR-Wp).str(int(prec*0.28),radius=False))
    E=mp.eigsy(Amp,eigvals_only=True)
    ev=sorted(float(e) for e in E)  # ALGEBRAIC order (signed)
    return ev
for c,N in [(13,10),(30,12)]:
    evf=build(c,N,1200,True); eva=build(c,N,1200,False)
    print(f"c={c} N={N}:")
    print(f"  FULL tau:   algebraic min-eig={evf[0]:+.4e}  (n_neg={sum(1 for e in evf if e<0)})  max={evf[-1]:.3f}")
    print(f"  arch+pole:  algebraic min-eig={eva[0]:+.4e}  (n_neg={sum(1 for e in eva if e<0)})  max={eva[-1]:.3f}")
    print(f"  -> arch+pole is {'INDEFINITE (has neg eigs)' if eva[0]<-1e-9 else 'PD'};  full is {'PD' if evf[0]>-1e-9 else 'INDEFINITE'}")

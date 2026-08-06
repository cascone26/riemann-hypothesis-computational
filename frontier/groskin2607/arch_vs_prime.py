"""Is the Weil-form saturation ARCHIMEDEAN (Connes/Sonin/prolate, understood) or PRIME-driven
(where RH lives, potentially our territory)? Decompose Groskin's tau = W02(pole) - WR(arch) - Wp(prime).
Compare min-eig of: arch+pole only (Wp=0) vs FULL. If arch-only already plunges, saturation is
archimedean (Connes). If the prime subtraction creates/drives the near-null, it's a prime phenomenon."""
import sys, math, numpy as np
sys.path.insert(0,'.')
from arb_ldlt_certify import (arb_closed_forms, arb_J, arb_kappa, prime_powers_up_to)
from flint import arb, acb, arb_mat, ctx
import mpmath as mp
def build(c,N,prec,include_prime=True):
    ctx.prec=prec
    S,CC,XC,L=arb_closed_forms(N,c,prec); PI=arb.pi(); sp2=16*PI*PI; l2=L*L
    pref02=32*L*(L/4).sinh()**2; kappa=arb_kappa(L); J=arb_J(L)
    pdata=prime_powers_up_to(c)
    weights=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    positions=[arb(q).log() for (q,p) in pdata]; M=len(weights)
    def Ss(nn): return S[nn] if nn>=0 else -S[-nn]
    DIM=2*N+1; A=arb_mat(DIM,DIM)
    for i in range(DIM):
        n=i-N
        for j in range(i,DIM):
            m=j-N; num=l2-sp2*m*n; den=(l2+sp2*m*m)*(l2+sp2*n*n); W02=pref02*num/den
            if n==m: WR=kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]
            else: WR=(Ss(m)-Ss(n))/(PI*(n-m))
            Wp=arb(0)
            if include_prime:
                for idx in range(M):
                    y=positions[idx]
                    if n==m: q=2*(1-y/L)*(2*PI*n*y/L).cos()
                    else: q=((2*PI*m*y/L).sin()-(2*PI*n*y/L).sin())/(PI*(n-m))
                    Wp+=weights[idx]*q
            val=W02-WR-Wp; A[i,j]=val; A[j,i]=val
    return A,DIM
def mineig(A,DIM,prec):
    mp.mp.dps=int(prec*0.29); Amp=mp.matrix(DIM,DIM)
    for i in range(DIM):
        for j in range(DIM): Amp[i,j]=mp.mpf(A[i,j].str(int(prec*0.28),radius=False))
    return min(float(mp.log10(abs(e))) for e in mp.eigsy(Amp,eigvals_only=True))
print(f"{'c':>5}{'N':>4}   {'min-eig ARCH+POLE only':>24}   {'min-eig FULL(+primes)':>22}   driver")
for c,N in [(13,12),(30,12),(30,16),(100,12),(100,16)]:
    prec=max(900,int((2.6*N+120)*4.2))
    Aa,D=build(c,N,prec,include_prime=False); va=mineig(Aa,D,prec)
    Af,_=build(c,N,prec,include_prime=True);  vf=mineig(Af,D,prec)
    driver="PRIME drives plunge" if vf<va-2 else "ARCHIMEDEAN (prolate/Connes)" if abs(vf-va)<2 else "both"
    print(f"{c:>5}{N:>4}   {va:>24.2f}   {vf:>22.2f}   {driver}  (Δ={vf-va:+.1f})")

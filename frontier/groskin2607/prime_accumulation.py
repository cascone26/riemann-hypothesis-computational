"""How many primes drive the saturation plunge? Compute min-eig of tau with the prime sum
restricted to {p<=P} for growing P. If p=2 alone (or few small primes) reproduces most of the
plunge, the RH-hard core is essentially LOW-DIMENSIONAL (a small-prime statement)."""
import sys, math, numpy as np
sys.path.insert(0,'.')
from arb_ldlt_certify import arb_closed_forms, arb_J, arb_kappa, prime_powers_up_to
from flint import arb, ctx
import mpmath as mp
def mineig_primesubset(c,N,prec,pmax_incl):
    ctx.prec=prec; S,CC,XC,L=arb_closed_forms(N,c,prec); PI=arb.pi(); sp2=16*PI*PI; l2=L*L
    pref02=32*L*(L/4).sinh()**2; kappa=arb_kappa(L); J=arb_J(L)
    pdata=[(q,p) for (q,p) in prime_powers_up_to(c) if p<=pmax_incl]
    weights=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]; positions=[arb(q).log() for (q,p) in pdata]; M=len(weights)
    def Ss(nn): return S[nn] if nn>=0 else -S[-nn]
    DIM=2*N+1; mp.mp.dps=int(prec*0.29); Amp=mp.matrix(DIM,DIM)
    for i in range(DIM):
        n=i-N
        for j in range(DIM):
            m=j-N; num=l2-sp2*m*n; den=(l2+sp2*m*m)*(l2+sp2*n*n); W02=pref02*num/den
            if n==m: WR=kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]
            else: WR=(Ss(m)-Ss(n))/(PI*(n-m))
            Wp=arb(0)
            for idx in range(M):
                y=positions[idx]
                if n==m: qq=2*(1-y/L)*(2*PI*n*y/L).cos()
                else: qq=((2*PI*m*y/L).sin()-(2*PI*n*y/L).sin())/(PI*(n-m))
                Wp+=weights[idx]*qq
            Amp[i,j]=mp.mpf((W02-WR-Wp).str(int(prec*0.28),radius=False))
    return min(float(mp.log10(abs(e))) for e in mp.eigsy(Amp,eigvals_only=True))
c,N=30,16; prec=1400
print(f"c={c} N={N}: min-eig as primes are added (log10):")
print(f"  no primes (arch only): {mineig_primesubset(c,N,prec,1):.2f}")
for P in [2,3,5,7,11,13,29]:
    v=mineig_primesubset(c,N,prec,P)
    print(f"  primes p<={P:2}: {v:.2f}")

"""For the extremal (min-eig) test function v of Groskin's tau(c,N), compute EACH prime's
contribution to <v, Wp v> = the prime term. Look for arithmetic structure (mu? sign pattern?
magnitude decay?) connecting the Weil extremal to our sign-theorem/rigidity prime-zero dictionary."""
import sys, math, numpy as np
sys.path.insert(0,'.')
from arb_ldlt_certify import arb_closed_forms, arb_J, arb_kappa, prime_powers_up_to
from flint import arb, arb_mat, ctx
import mpmath as mp
from sympy import mobius, factorint
def extremal_v(c,N,prec):
    # build full tau, get min-eig eigenvector (float via mpmath)
    ctx.prec=prec; S,CC,XC,L=arb_closed_forms(N,c,prec); PI=arb.pi(); sp2=16*PI*PI; l2=L*L
    pref02=32*L*(L/4).sinh()**2; kappa=arb_kappa(L); J=arb_J(L)
    pdata=prime_powers_up_to(c)
    weights=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]; positions=[arb(q).log() for (q,p) in pdata]; M=len(weights)
    def Ss(nn): return S[nn] if nn>=0 else -S[-nn]
    DIM=2*N+1
    mp.mp.dps=int(prec*0.29); Amp=mp.matrix(DIM,DIM)
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
    E,V=mp.eigsy(Amp); idx0=min(range(DIM),key=lambda k:abs(float(E[k])))
    v=[mp.mpf(V[i,idx0]) for i in range(DIM)]
    # per-prime contribution to <v,Wp v>
    Lf=float(mp.log(c)); PIf=math.pi
    print(f"c={c} N={N}  min-eig=10^{float(mp.log10(abs(E[idx0]))):.1f}")
    print(f"  {'q=p^a':>7} {'Lambda/sqrt(q)':>14} {'<v,Wp_q v>':>13} {'sign':>5} {'mu(q)':>6}")
    contribs=[]
    for idx,(q,p) in enumerate(pdata):
        y=float(positions[idx]); w=float(weights[idx]); s=mp.mpf(0)
        for i in range(DIM):
            n=i-N
            for j in range(DIM):
                m=j-N
                if n==m: qq=2*(1-y/Lf)*math.cos(2*PIf*n*y/Lf)
                else: qq=(math.sin(2*PIf*m*y/Lf)-math.sin(2*PIf*n*y/Lf))/(PIf*(n-m))
                s+=v[i]*v[j]*qq
        contrib=w*float(s); mu=int(mobius(q)) if q==p else 0  # mu only for primes; prime powers mu=0
        contribs.append((q,p,contrib))
        print(f"  {q:>7} {w:>14.5f} {contrib:>13.3e} {'+' if contrib>0 else '-':>5} {mu:>6}")
    tot=sum(c[2] for c in contribs); print(f"  TOTAL <v,Wp v> = {tot:.4e}  (= prime term; should ~ arch term for extremal)")
extremal_v(30,12,1200)

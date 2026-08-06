"""CORRECT structural decomposition of Groskin's tau = W02(pole) - WR(arch) - Wp(prime).
Algebraic (signed) min-eig of every sub-combination, to see which pieces are pos/indef and
how positivity actually arises. Fixes the §741 min-|eig| bug."""
import sys, numpy as np
sys.path.insert(0,'.')
from arb_ldlt_certify import arb_closed_forms, arb_J, arb_kappa, prime_powers_up_to
from flint import arb, ctx
import mpmath as mp
def mats(c,N,prec):
    ctx.prec=prec; S,CC,XC,L=arb_closed_forms(N,c,prec); PI=arb.pi(); sp2=16*PI*PI; l2=L*L
    pref02=32*L*(L/4).sinh()**2; kappa=arb_kappa(L); J=arb_J(L)
    pdata=prime_powers_up_to(c); weights=[arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    positions=[arb(q).log() for (q,p) in pdata]; M=len(weights)
    def Ss(nn): return S[nn] if nn>=0 else -S[-nn]
    DIM=2*N+1; W02=np.zeros((DIM,DIM)); WR=np.zeros((DIM,DIM)); Wp=np.zeros((DIM,DIM))
    for i in range(DIM):
        n=i-N
        for j in range(DIM):
            m=j-N; num=l2-sp2*m*n; den=(l2+sp2*m*m)*(l2+sp2*n*n)
            W02[i,j]=float((pref02*num/den).mid() if hasattr(pref02*num/den,'mid') else pref02*num/den)
            wr=(kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]) if n==m else (Ss(m)-Ss(n))/(PI*(n-m))
            WR[i,j]=float(wr.mid() if hasattr(wr,'mid') else wr)
            wp=arb(0)
            for idx in range(M):
                y=positions[idx]
                if n==m: qq=2*(1-y/L)*(2*PI*n*y/L).cos()
                else: qq=((2*PI*m*y/L).sin()-(2*PI*n*y/L).sin())/(PI*(n-m))
                wp+=weights[idx]*qq
            Wp[i,j]=float(wp.mid() if hasattr(wp,'mid') else wp)
    return W02,WR,Wp
def me(M): 
    e=np.linalg.eigvalsh(M); return e[0],e[-1]  # algebraic min,max
c,N=30,12
P,A,Q=mats(c,N,900)  # P=pole W02, A=arch WR, Q=prime Wp
combos={"pole W02":P, "arch WR":A, "prime Wp":Q, "pole-arch (W02-WR)":P-A,
        "pole-prime (W02-Wp)":P-Q, "-(arch)-(prime)=-(WR+Wp)":-(A+Q),
        "arch+prime (WR+Wp)":A+Q, "FULL (W02-WR-Wp)":P-A-Q}
print(f"c={c} N={N}: algebraic [min-eig, max-eig] of each piece/combo:")
for name,M in combos.items():
    mn,mx=me(M); tag="PD" if mn>-1e-9 else ("NEG-def" if mx<1e-9 else "INDEF")
    print(f"  {name:28s} [{mn:+.4e}, {mx:+.4e}]  {tag}")
print("\ncommutators (normalized ||[X,Y]||/||X||||Y||):")
for (n1,X),(n2,Y) in [(("pole",P),("arch",A)),(("pole",P),("prime",Q)),(("arch",A),("prime",Q))]:
    cn=np.linalg.norm(X@Y-Y@X)/(np.linalg.norm(X)*np.linalg.norm(Y))
    print(f"  [{n1},{n2}] = {cn:.4f}")

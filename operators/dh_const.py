import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
import mpmath as mp
mp.mp.dps=18
from flint import arb, ctx
import dh_certified_tau as DH   # scaffold closed_forms_nu (exact Arb)

# numeric psi_R'(n) at nu=3/4, cond=5 (validated machinery)
def psiR(x,L,nu,cond):
    x=mp.mpf(x);nu=mp.mpf(nu);L=mp.mpf(L);logc=mp.log(mp.pi/cond)
    hp=lambda r: mp.re(mp.digamma(nu+1j*r/2))-logc
    def Sr(r):
        B=2*mp.pi*x/L;den=B*B-r*r
        if abs(den)<mp.mpf(10)**-18: r=r+mp.mpf(10)**-12;den=B*B-r*r
        return (mp.cos(r*L)-mp.cos(2*mp.pi*x))*B/den
    f=lambda r: hp(r)*Sr(r);B=abs(2*mp.pi*x/L)
    per=mp.pi/L;R=mp.mpf(40);pts=[mp.mpf(0)];k=2
    while k*per<R: pts.append(k*per);k+=2
    pts.append(R);allp=sorted(set([float(p) for p in pts]+[float(B)]))
    val=2*mp.quad(f,allp)+2*B*mp.cos(2*mp.pi*x)*mp.quad(lambda r: hp(r)/r**2,[R,4*R,mp.inf])
    return -val/(2*mp.pi**2)

c=14;N=4;prec=400;ctx.prec=prec
nu=mp.mpf(3)/4
S,CC,XC,L=DH.closed_forms_nu(N,c,prec,0.75)   # exact-Arb scaffold at nu=3/4
Lf=float(L.mid()); Lm=float(mp.log(c))
h=mp.mpf(10)**-5
# numeric diag = (kappa+J) + 2 CC(n) - (2/L) XC(n)  =>  (kappa+J) = diag - 2 CC(n) + (2/L) XC(n)
print("(kappa+J)_DH from numeric diag + exact-Arb scaffold CC/XC at nu=3/4, cond=5:",flush=True)
for n in [1,2,3]:
    diag=float((psiR(n+h,Lm,nu,5)-psiR(n-h,Lm,nu,5))/(2*h))
    ccn=float(CC[abs(n)].mid()); xcn=float(XC[abs(n)].mid())
    KJ=diag - 2*ccn + (2/Lf)*xcn
    print(f"  n={n}: numeric_diag={diag:+.6e}  CC={ccn:+.5e} XC={xcn:+.5e}  (kappa+J)_DH={KJ:+.6e}",flush=True)
# zeta cross-check: same formula with nu=1/4,cond=1 must recover Groskin's exact kappa+J
import arb_ldlt_certify as G
kapz=float(G.arb_kappa(L).mid()); Jz=float(G.arb_J(L).mid())
print(f"\nZETA cross-check: exact Groskin kappa+J = {kapz+Jz:+.6e}",flush=True)
Sz,CCz,XCz,Lz=G.arb_closed_forms(N,c,prec)
for n in [1,2,3]:
    diag=float((psiR(n+h,Lm,mp.mpf(1)/4,1)-psiR(n-h,Lm,mp.mpf(1)/4,1))/(2*h))
    KJ=diag - 2*float(CCz[abs(n)].mid()) + (2/float(Lz.mid()))*float(XCz[abs(n)].mid())
    print(f"  n={n}: numeric-derived (kappa+J)_zeta={KJ:+.6e}  (target {kapz+Jz:+.6e}, err {abs(KJ-(kapz+Jz)):.1e})",flush=True)

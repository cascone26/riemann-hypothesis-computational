import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
import mpmath as mp
mp.mp.dps=15
from flint import arb, ctx
import arb_ldlt_certify as G

def psiR_fast(x, L, nu, cond):
    x=mp.mpf(x); nu=mp.mpf(nu); L=mp.mpf(L)
    logc=mp.log(mp.pi/cond)
    hplus=lambda r: mp.re(mp.digamma(nu+1j*r/2))-logc
    def Sr(r):
        B=2*mp.pi*x/L; den=B*B-r*r
        if abs(den)<mp.mpf(10)**-18: r=r+mp.mpf(10)**-12; den=B*B-r*r
        return (mp.cos(r*L)-mp.cos(2*mp.pi*x))*B/den
    f=lambda r: hplus(r)*Sr(r)
    B=abs(2*mp.pi*x/L) if x!=0 else mp.mpf(1)
    per=mp.pi/L; R=mp.mpf(40)
    pts=[mp.mpf(0)]; k=4
    while k*per<R: pts.append(k*per); k+=4
    pts.append(R)
    allp=sorted(set([float(p) for p in pts]+[float(B)]))
    val=2*mp.quad(f,allp)
    tail=mp.quad(lambda r: hplus(r)/r**2,[R,4*R,mp.inf])
    val+=2*B*mp.cos(2*mp.pi*x)*tail
    return -val/(2*mp.pi**2)

c=14;N=6;L=float(mp.log(c)); ctx.prec=300
S,CC,XC,La=G.arb_closed_forms(N,c,300); kap=G.arb_kappa(La); J=G.arb_J(La); PI=arb.pi()
def clo(n,m):
    Ss=lambda k: S[k] if k>=0 else -S[-k]
    if n==m: return float((kap+2*CC[abs(n)]+J-(2/La)*XC[abs(n)]).mid())
    return float(((Ss(m)-Ss(n))/(PI*(n-m))).mid())
h=mp.mpf(10)**-5
print("DIAG (tail-corrected,dps15,coarse):",flush=True)
for n in [1,2,3]:
    d=float((psiR_fast(n+h,L,mp.mpf(1)/4,1)-psiR_fast(n-h,L,mp.mpf(1)/4,1))/(2*h))
    print(f"  n={n}: num={d:+.5e} clo={clo(n,n):+.5e} relerr={abs(d-clo(n,n))/abs(clo(n,n)):.1e}",flush=True)
v={k:psiR_fast(k,L,mp.mpf(1)/4,1) for k in range(4)}
print("OFFDIAG:",flush=True)
for (m,n) in [(1,0),(2,1),(3,1)]:
    nu2=float((v[m]-v[n])/(m-n)); print(f"  (m{m},n{n}): num={nu2:+.5e} clo={clo(n,m):+.5e} relerr={abs(nu2-clo(n,m))/(abs(clo(n,m))+1e-30):.1e}",flush=True)

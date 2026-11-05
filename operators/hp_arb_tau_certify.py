"""S2880: RIGOROUS ball-arithmetic build of tau = W02 - W_R - W_p (CCM = Shi S_NL) and LDL^T certificate of tau - e*I > 0.
python-flint (Arb). Arch integrals: acb.integral on [d, L] + explicit bound on [0, d]."""
import sys
from flint import arb, acb, arb_mat, ctx
lam2s=sys.argv[1]; N=int(sys.argv[2]); prec=int(sys.argv[3]); eshift=sys.argv[4]
ctx.prec=prec
from fractions import Fraction as _F; _q=_F(lam2s); LAM2=arb(_q.numerator)/arb(_q.denominator); L=LAM2.log(); lam2=lam2s; PI=arb.pi(); d=arb(10)**(-60)
M=2*N+1; idx=list(range(-N,N+1))
# von Mangoldt k < lam2 (k = lam2 has weight factor (L - log k)=0 / sin terms vanish anyway)
def vm(k):
    for p in range(2,k+1):
        if k%p==0:
            kk=k
            while kk%p==0: kk//=p
            return arb(p).log() if kk==1 else None
pw=[(k,vm(k)) for k in range(2,int(_F(lam2s))+1) if vm(k) is not None and k<_F(lam2s)]
def Jsin(n):
    if n==0: return arb(0)
    a=2*PI*n/L
    f=lambda x,an: acb(x/2).exp()*(a*x).sin()/(acb(x).exp()-acb(-x).exp())
    I=acb.integral(f,acb(d),acb(L),rel_tol=arb(2)**(-prec+20),eval_limit=10**7,deg_limit=prec).real
    # [0,d]: |e^{x/2} sin(ax)/(2 sinh x)| <= a*e^{d}/2
    return I+arb(0,(a*d.exp()*d/2).upper())
def Ddiag(n):
    a=2*PI*n/L
    f=lambda x,an: (acb(x/2).exp()*2*(L-x)*(a*x).cos()/L-2)/(acb(x).exp()-acb(-x).exp())
    I=acb.integral(f,acb(d),acb(L),rel_tol=arb(2)**(-prec+20),eval_limit=10**7,deg_limit=prec).real
    K=(1+2/L+a*a*L+2*a)*arb(1).exp()*4   # crude sup bound of integrand on [0,d]
    return I+arb(0,(K*d).upper())
Js={n:Jsin(n) for n in range(N+1)}
print('max rad J:',max(float(Js[n].rad()) for n in Js),flush=True)
Jsig=lambda n: Js[abs(n)]*(1 if n>=0 else -1)
Dd={n:Ddiag(n) for n in range(N+1)}
print('max rad D:',max(float(Dd[n].rad()) for n in Dd),flush=True)
archconst=arb.const_euler()+(4*PI*(L.exp()-1)/(L.exp()+1)).log()
def omega(n,m,y):
    if n==m: return 2*(L-y)*(2*PI*n*y/L).cos()/L
    return ((2*PI*m*y/L).sin()-(2*PI*n*y/L).sin())/(PI*(n-m))
sh2=(L/4).sinh()**2
T=arb_mat(M,M)
for i,n in enumerate(idx):
    for j,m in enumerate(idx):
        W02=32*L*sh2*(L*L-16*PI*PI*m*n)/((L*L+16*PI*PI*m*m)*(L*L+16*PI*PI*n*n))
        WR=(archconst+Dd[abs(n)]) if n==m else (Jsig(m)-Jsig(n))/(PI*(n-m))
        WP=sum((lk*arb(k)**arb(-0.5)*omega(n,m,arb(k).log()) for k,lk in pw), arb(0))
        T[i,j]=W02-WR-WP-(arb(eshift) if i==j else arb(0))
# RIGOROUS certificate: exact-rational LDL^T on midpoints + Weyl bound for the radii.
from fractions import Fraction
DIG=prec*3//10
def tofrac(x):
    st=x.mid().str(DIG,radius=False)
    return Fraction(st.replace(' ','')) if 'e' not in st else Fraction(st)
maxrad=max(float(T[i,j].rad()) for i in range(M) for j in range(M))
rho=Fraction(M)*Fraction(maxrad)*2 + Fraction(1,10**(DIG-5))*M  # radius (Weyl, ||.||_2<=M*max) + decimal rounding
A=[[tofrac(T[i,j]) for j in range(M)] for i in range(M)]
for i in range(M): A[i][i]-=rho
ok=True; pivots=[]
for k in range(M):
    p=A[k][k]; pivots.append(p)
    if p<=0: ok=False; break
    for i2 in range(k+1,M):
        f=A[i2][k]/p
        if f:
            for j in range(k+1,M): A[i2][j]-=f*A[k][j]
mn=min(pivots)
print(f"lam2={lam2} N={N} prec={prec} shift e={eshift}: CERTIFIED (tau - eI) >= (mid - rho - e) > 0 : {ok}  pivots ok={len(pivots) if ok else len(pivots)-1}/{M}  min pivot~{float(mn):.4e}  rho={float(rho):.2e}")

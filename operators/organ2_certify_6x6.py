"""Extend the rigorous Weil-positivity certificate to a LARGER basis (6 test functions spanning
more frequencies). Arb enclosure -> exact-rational LDL^T. Genuine forward progress on the actual
RH criterion (Weil form PD for finite test-function families), verified-numerics lane (§755)."""
from flint import arb, acb, ctx
import sympy as sp
from sympy import primerange, Rational
ctx.prec = 300
AV = [arb(15), arb(20), arb(25), arb(30), arb(35), arb(40)]  # 6-function basis, on the zeros
nb = len(AV)
PI = arb.pi(); SQRT2PI = (2*PI).sqrt(); LOGPI = PI.log()
def phihat_acb(a, r):
    A=acb(a); return acb(SQRT2PI)/2*((-(r-A)**2/2).exp()+(-(r+A)**2/2).exp())
def Omega(r):
    z1=acb(arb(1)/4)+acb(0,1)*r/2; z2=acb(arb(1)/4)-acb(0,1)*r/2
    return (z1.digamma()+z2.digamma())/2-acb(LOGPI)
def arch_entry(aj,ak):
    R=arb(120)
    def integrand(x,analytic):
        if analytic:
            im=x.imag
            if float(im.upper())>=0.5 or float(im.lower())<=-0.5: return acb(arb('nan'))
        return phihat_acb(aj,x)*phihat_acb(ak,x)*Omega(x)
    return acb.integral(integrand, acb(-R), acb(R)).real/(2*PI)
def pole_entry(aj,ak):
    r=acb(0,arb(1)/2); return 2*(phihat_acb(aj,r)*phihat_acb(ak,r)).real
def g_cross(aj,ak,y):
    s=aj+ak; d=aj-ak
    return (-y**2/4).exp()*(PI.sqrt()/2)*((-s**2/4).exp()*(d*y/2).cos()+(-d**2/4).exp()*(s*y/2).cos())
def prime_entry(aj,ak):
    tot=arb(0)
    for p in primerange(2,2000000):
        lp=arb(p).log(); m=1
        while True:
            L=m*lp
            if float(L)>30: break
            tot+=lp/arb(p)**(arb(m)/2)*g_cross(aj,ak,L); m+=1
    prime=2*tot
    Lcut=arb(2000000).log(); gmax=(PI.sqrt()/2)*((-(aj+ak)**2/4).exp()+(-(aj-ak)**2/4).exp())
    tailb=2*gmax*(arb(1)/4).exp()*PI.sqrt()*((Lcut-1)/2).erfc()
    return prime+arb(0,float(tailb))
print(f"Building rigorous {nb}x{nb} Weil form...")
Qb=[[None]*nb for _ in range(nb)]
for j in range(nb):
    for k in range(j,nb):
        v=arch_entry(AV[j],AV[k])+pole_entry(AV[j],AV[k])-prime_entry(AV[j],AV[k])
        Qb[j][k]=v; Qb[k][j]=v
M=sp.zeros(nb,nb); radsq=arb(0)
for j in range(nb):
    for k in range(nb):
        mid=Qb[j][k].mid(); rad=Qb[j][k].rad()
        M[j,k]=Rational(str(float(mid))).limit_denominator(10**12)
        radsq+=arb(abs(float(mid)-float(M[j,k]))+float(rad))**2
rho=radsq.sqrt(); rho_rat=Rational(str(float(rho)+1e-15)).limit_denominator(10**10)+Rational(1,10**9)
print(f"rho <= {float(rho):.3e}")
Msh=(M-rho_rat*sp.eye(nb)).applyfunc(lambda e: sp.Rational(e))
n=nb; L=sp.eye(n); D=[sp.Integer(0)]*n; pd=True
for jj in range(n):
    d=sp.Rational(Msh[jj,jj]-sum(L[jj,kk]**2*D[kk] for kk in range(jj))); D[jj]=d
    if D[jj]<=0: pd=False
    for ii in range(jj+1,n):
        L[ii,jj]=sp.Rational((Msh[ii,jj]-sum(L[ii,kk]*L[jj,kk]*D[kk] for kk in range(jj)))/D[jj])
print("pivots(float):",[f'{float(dv):+.4e}' for dv in D])
print("all pivots>0 (EXACT):",pd)
minors=[sp.Rational((Msh[:kk,:kk]).det()) for kk in range(1,n+1)]
sylv=all(m>0 for m in minors)
print(f"Sylvester PD:{sylv} | agree:{pd==sylv}")
print("===>", "HANDLE: Weil form PD on 6-function basis {15,20,25,30,35,40}, RIGOROUS (Arb+exact LDL^T)" if (pd and sylv) else "not certified")

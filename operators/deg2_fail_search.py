# Degree-2 FAIL control search: weight-24 non-eigenform combinations g = Delta*E12 + lam*Delta^2.
# Scan for LOW off-line zeros (argument principle total > on-line count) = a real degree-2 D-H analog.
import mpmath as mp
from fractions import Fraction as Fr
import numpy as np
mp.mp.dps=25
NM=60
# tau
P=[0]*(NM+1);P[0]=1
for n in range(1,NM+1):
    for k in range(NM,n-1,-1):P[k]-=P[k-n]
def pmul(a,b):
    c=[0]*(NM+1)
    for k in range(NM+1):
        if a[k]:
            for j in range(NM+1-k):c[k+j]+=a[k]*b[j]
    return c
Qp=[0]*(NM+1);Qp[0]=1;base=P[:];e=24
while e:
    if e&1:Qp=pmul(Qp,base)
    e>>=1
    if e:base=pmul(base,base)
tau=[0]+[Qp[n-1] for n in range(1,NM+1)]
# Delta^2 coeffs: a_n(D^2)=sum tau(i)tau(n-i)
d2=[0]*(NM+1)
for n in range(2,NM+1):
    d2[n]=sum(tau[i]*tau[n-i] for i in range(1,n))
# Delta*E12
def sig11(n):
    s=0;dd=1
    while dd*dd<=n:
        if n%dd==0:
            s+=dd**11
            if dd!=n//dd:s+=(n//dd)**11
        dd+=1
    return s
c12=[Fr(0)]*(NM+1);c12[0]=Fr(1)
for n in range(1,NM+1):c12[n]=Fr(65520,691)*sig11(n)
de=[Fr(0)]*(NM+1)
for n in range(1,NM+1):
    de[n]=sum(Fr(tau[i])*c12[n-i] for i in range(1,n+1))
K=24;TP=2*mp.pi
def scan(lam):
    a=[mp.mpf(0)]+[mp.mpf(de[n].numerator)/de[n].denominator + lam*d2[n] for n in range(1,NM+1)]
    def Lam(s):
        s=mp.mpc(s);t=mp.mpc(0)
        for n in range(1,NM+1):
            x=TP*n;t+=a[n]*((x)**(-s)*mp.gammainc(s,x)+(x)**(s-K)*mp.gammainc(K-s,x))
        return t
    # on-line count [0,36]
    ts=np.arange(0.3,36,0.2);prev=None;onc=0
    for tt in ts:
        v=float(mp.re(Lam(mp.mpc(12,tt))))
        if prev is not None and prev*v<0:onc+=1
        prev=v
    # argument principle total, strip [10.5,13.5]x[0,36]
    def cnt(t0,t1,lo=10.5,hi=13.5,ns=500):
        T=np.linspace(t0,t1,ns);S=np.linspace(hi,lo,ns)
        pts=[mp.mpc(hi,t) for t in T]+[mp.mpc(s,t1) for s in S]+[mp.mpc(lo,t) for t in T[::-1]]+[mp.mpc(s,t0) for s in S[::-1]]
        tot=mp.mpf(0);vp=Lam(pts[0])
        for p in pts[1:]:
            v=Lam(p);tot+=mp.arg(v/vp);vp=v
        return float(tot/(2*mp.pi))
    tot=cnt(0.3,36)
    return onc,tot
for lam in [0.0, 0.5, 2.0, -1.0]:
    onc,tot=scan(lam)
    off=tot-onc
    print(f"lam={lam:+.1f}: on-line={onc}  total(arg-principle)={tot:.2f}  off-line~{off:.2f}  {'<-- OFF-LINE ZEROS!' if off>0.5 else 'all on-line'}",flush=True)

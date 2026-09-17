import mpmath as mp
from fractions import Fraction as Fr
import numpy as np
mp.mp.dps=52
NM=32
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
def sig11(n):
    s=0;d=1
    while d*d<=n:
        if n%d==0:
            s+=d**11
            if d!=n//d:s+=(n//d)**11
        d+=1
    return s
c12=[Fr(0)]*(NM+1);c12[0]=Fr(1)
for n in range(1,NM+1):c12[n]=Fr(65520,691)*sig11(n)
a=[Fr(0)]*(NM+1)
for n in range(1,NM+1):
    s=Fr(0)
    for i in range(1,n+1):s+=tau[i]*c12[n-i]
    a[n]=s
av=[mp.mpf(0)]+[mp.mpf(a[n].numerator)/a[n].denominator for n in range(1,NM+1)]
K=24;TP=2*mp.pi
def ReLam(t):
    s=mp.mpc(12,t);tot=mp.mpc(0)
    for n in range(1,NM+1):
        x=TP*n;tot+=av[n]*((x)**(-s)*mp.gammainc(s,x)+(x)**(s-K)*mp.gammainc(K-s,x))
    return float(mp.re(tot))
# compare grids on [80,90] (smaller window, faster): does finer grid recover zeros?
for dt in [0.05,0.02,0.01]:
    ts=np.arange(80.0,90.0,dt);prev=None;c=0
    for t in ts:
        v=ReLam(t)
        if prev is not None and prev*v<0:c+=1
        prev=v
    print(f"[80,90] grid={dt}: on-line sign-changes = {c}",flush=True)

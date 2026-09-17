# Exact completed L-function of Delta*E12 (weight 24, eps=+1) via incomplete gamma; locate off-line zeros.
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps=30
NM=80
# tau(n) exact
P=[0]*(NM+1); P[0]=1
for n in range(1,NM+1):
    for k in range(NM,n-1,-1): P[k]-=P[k-n]
def pmul(a,b):
    c=[0]*(NM+1)
    for k in range(NM+1):
        if a[k]:
            for j in range(NM+1-k): c[k+j]+=a[k]*b[j]
    return c
Qp=[0]*(NM+1);Qp[0]=1;base=P[:];e=24
while e:
    if e&1:Qp=pmul(Qp,base)
    e>>=1
    if e:base=pmul(base,base)
tau=[0]*(NM+1)
for n in range(1,NM+1): tau[n]=Qp[n-1]
def sig11(n):
    s=0;d=1
    while d*d<=n:
        if n%d==0:
            s+=d**11
            if d!=n//d: s+=(n//d)**11
        d+=1
    return s
c12=[Fr(0)]*(NM+1);c12[0]=Fr(1)
for n in range(1,NM+1): c12[n]=Fr(65520,691)*sig11(n)
a=[Fr(0)]*(NM+1)
for n in range(1,NM+1):
    s=Fr(0)
    for i in range(1,n+1): s+=tau[i]*c12[n-i]
    a[n]=s
av=[mp.mpf(0)]+[mp.mpf(a[n].numerator)/a[n].denominator for n in range(1,NM+1)]
K=24; TWO_PI=2*mp.pi
def Lam(s):   # completed L, weight 24, eps=+1; real on Re s=12
    s=mp.mpc(s); tot=mp.mpc(0)
    for n in range(1,NM+1):
        x=TWO_PI*n
        tot+=av[n]*( (x)**(-s)*mp.gammainc(s,x) + (x)**(s-K)*mp.gammainc(K-s,x) )
    return tot
# 1) on-line zeros: Lam(12+it) real; find sign changes
print("on-line scan Lam(12+it), t in [0,40]:",flush=True)
import numpy as np
ts=np.arange(0.2,40,0.2); prev=None; onzeros=[]
for t in ts:
    v=mp.re(Lam(12+1j*t))
    if prev is not None and prev*v<0:
        z=mp.findroot(lambda tt: mp.re(Lam(12+1j*tt)), t-0.1); onzeros.append(float(z))
    prev=v
print(f"  {len(onzeros)} on-line zeros up to t=40: {[round(z,3) for z in onzeros[:20]]}",flush=True)
# 2) argument principle: count ALL zeros in box Re in [11,13], t in [0,T] ; compare to on-line count
def count_box(t0,t1,sig_lo=11.0,sig_hi=13.0,ns=400):
    # integrate d(arg Lam) around rectangle
    pts=[]
    T=np.linspace(t0,t1,ns)
    for t in T: pts.append(mp.mpc(sig_hi,t))
    S=np.linspace(sig_hi,sig_lo,ns)
    for s in S: pts.append(mp.mpc(s,t1))
    for t in T[::-1]: pts.append(mp.mpc(sig_lo,t))
    for s in S[::-1]: pts.append(mp.mpc(s,t0))
    tot=mp.mpf(0); vprev=Lam(pts[0])
    for p in pts[1:]:
        v=Lam(p); tot+=mp.arg(v/vprev); vprev=v
    return tot/(2*mp.pi)
for (t0,t1) in [(0.1,10),(10,20),(20,30),(30,40)]:
    nbox=count_box(t0,t1); non=sum(1 for z in onzeros if t0<z<=t1)
    print(f"  box t in [{t0},{t1}]: total zeros ~ {float(nbox):+.2f}, on-line = {non}  => off-line ~ {float(nbox)-non:.2f}",flush=True)

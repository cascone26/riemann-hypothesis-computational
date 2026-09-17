# Direct off-line zero search: scan |Lambda(sigma+it)| off the critical line Re=12 in [80,90].
# A near-zero minimum at sigma != 12 = DIRECT proof of an off-line zero (unlike a count deficit).
import mpmath as mp
from fractions import Fraction as Fr
import numpy as np
mp.mp.dps=48
NM=30
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
def Lam(s):
    s=mp.mpc(s);tot=mp.mpc(0)
    for n in range(1,NM+1):
        x=TP*n;tot+=av[n]*((x)**(-s)*mp.gammainc(s,x)+(x)**(s-K)*mp.gammainc(K-s,x))
    return tot
# normalize by |Lambda(13.5+it)| scale? No -- compare log|Lam| across sigma at each t to spot off-line dips.
# For each t, find min over sigma of |Lam|; an off-line zero shows |Lam|->0 at sigma!=12.
ts=np.arange(80.0,90.01,0.1)
sigs=[12.0,12.4,12.8,13.2,13.6,14.0]
print(f"{'t':>6} " + " ".join(f"s={s:>5.2f}" for s in sigs) + "   argmin(off-line?)",flush=True)
offline_hits=[]
for t in ts:
    row=[float(mp.fabs(Lam(mp.mpc(s,t)))) for s in sigs]
    lr=[np.log10(r) if r>0 else -99 for r in row]
    # off-line dip: some sigma>12 has |Lam| comparable-or-smaller than on-line value AND is a local min
    j=int(np.argmin(row))
    flag=""
    if sigs[j]>12.01 and row[j]<0.2*row[0]:  # off-line value notably smaller than on-line
        flag=f"  <== off-line dip at s={sigs[j]:.2f} (|L|={row[j]:.1e} vs on-line {row[0]:.1e})"
        offline_hits.append((round(float(t),2),sigs[j],row[j],row[0]))
    print(f"{t:6.2f} " + " ".join(f"{v:6.1f}" for v in lr) + flag,flush=True)
print(f"\noff-line dips found: {len(offline_hits)}",flush=True)
for h in offline_hits[:10]: print("  ",h,flush=True)

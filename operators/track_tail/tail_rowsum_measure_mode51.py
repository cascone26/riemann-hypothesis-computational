import numpy as np
from scipy.special import digamma, zeta as hurwitz
import cmath, math

A = 0.549
b = math.pi/(2*A); log2 = math.log(2); twoa = 2*A; wp = log2/math.sqrt(2)
Om0 = digamma(0.25) - math.log(math.pi)
Phi0 = hurwitz(2, 0.25)
NL = 4000

def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3

def II(p,al,om,sc,t0,t1):
    s = complex(-al,om)
    F = lambda tt: cmath.exp(s*tt)*Pp(p,tt,s)
    v = F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag

def dterms(k,l):
    if k==l: return [(2*A*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
    al = k*l*b*(1.0/(k-l)-1.0/(k+l)); ak = -k*l*b*(1.0/(k-l)+1.0/(k+l))
    return [(al,0,l*b,'sin'),(ak,0,k*b,'sin')]

gt = [(-4,0,-0.5),(-4,0,0.5),(8,0,0.0),(-Om0/2,1,0.0),(-Phi0/4,0,0.0)]
for n in range(NL): gt.append((0.25/(n+0.25)**2,0,2*n+0.5))

def Qraw(k,l):
    tot = 0.0
    D = dterms(k,l)
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in D: tot += ag*ad*II(pg+pd,alp,om,sc,0,twoa)
    for (ad,pd,om,sc) in D:
        tot += wp*ad*II(pd+1,0,om,sc,log2,twoa)+(-wp*log2)*ad*II(pd,0,om,sc,log2,twoa)
    return tot/A

KMAX = 51
ODD = list(range(1,KMAX+1,2)); n=len(ODD)
Q = np.zeros((n,n))
for i,k in enumerate(ODD):
    for j,l in enumerate(ODD):
        if j<i: Q[i,j]=Q[j,i]
        else: Q[i,j]=Qraw(k,l)

# sanity vs NL=24000 table: diagonal + Gershgorin margins for k>23
print("sanity diag:", {k: round(Q[i,i],4) for i,k in enumerate(ODD) if k in (13,23,27)})
print()
print("row structure (fixed k, |Q[k,l]| vs l):")
for k in (15, 27, 39, 51):
    i = ODD.index(k)
    offs = [(l, abs(Q[i,ODD.index(l)])) for l in ODD if l!=k]
    rs = sum(v for _,v in offs)
    print(f"k={k:2d} diag={Q[i,i]:.4f} rowsum={rs:.4f} margin={Q[i,i]-rs:.4f}")
    big = sorted(offs, key=lambda t:-t[1])[:6]
    print("   largest:", " ".join(f"l={l}:{v:.2e}" for l,v in big))
# decay in |k-l| along a row
i = ODD.index(39)
print()
print("decay along row k=39: |Q[39,l]| by |k-l|:")
for l in ODD:
    if l!=39: print(f"  l={l:2d} |k-l|={abs(39-l):2d}  {abs(Q[i,ODD.index(l)]):.3e}")

print()
print("full Gershgorin sweep, window l<=51:")
mm=[]
for k in ODD:
    if k<15: continue
    i=ODD.index(k)
    rs=sum(abs(Q[i,ODD.index(l)]) for l in ODD if l!=k)
    mm.append((k,Q[i,i]-rs))
for k,g in mm: print(f"  k={k:2d} margin={g:.4f}")
print("min margin k in [15,51]:", min(g for _,g in mm))
np.save('/Users/scones/tmp/Q51.npy', Q)

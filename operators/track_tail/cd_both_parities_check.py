# C_d derivation check. Claim: R(k) := Q[k,k] - [log(kb/2pi) - amp*cos(kb*log2)] satisfies
#   R(k)*k = -1/2 - (2wp/pi)*sin(kb*log2) + O(1/k)     [both parities: Si sign flip is O(1/k^2)]
# so |R*k| <= 1/2 + 2wp/pi = 0.81205 + eps_k.  (ab = pi/2 exactly -> jump coeff wp/(ab) = 2wp/pi.)
import sys, math
import numpy as np
sys.path.insert(0,'/Users/scones/projects/riemann/rigor')
import mpmath as mp
mp.mp.dps = 25
A = mp.mpf('0.549'); b=mp.pi/(2*A); log2=mp.log(2); wp=log2/mp.sqrt(2)
Om0=mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi); Phi0=mp.zeta(2,mp.mpf(1)/4); NL=3000
def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3
def II(p,al,om,sc,t0,t1):
    s=mp.mpc(-al,om); F=lambda tt: mp.e**(s*tt)*Pp(p,tt,s); v=F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
def dterms_diag(k):
    return [(2*A*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
def Qdiag(k):
    tot=mp.mpf(0); D=dterms_diag(k)
    gt=[(-4,0,mp.mpf(-1)/2),(-4,0,mp.mpf(1)/2),(8,0,mp.mpf(0)),(-Om0/2,1,mp.mpf(0)),(-Phi0/4,0,mp.mpf(0))]
    for n in range(NL): gt.append((mp.mpf(1)/4/(n+mp.mpf(1)/4)**2,0,2*n+mp.mpf(1)/2))
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in D: tot+=ag*ad*II(pg+pd,alp,om,sc,0,2*A)
    for (ad,pd,om,sc) in D:
        tot+=wp*ad*II(pd+1,0,om,sc,log2,2*A)+(-wp*log2)*ad*II(pd,0,om,sc,log2,2*A)
    return float(tot/A)
amp=float(wp*(2*A-log2)/A); env=0.5+float(2*wp/mp.pi)
print(f"predicted envelope |R*k| <= 1/2 + 2wp/pi = {env:.5f}")
print("parity  k    Q[k,k]    R*k       model=-1/2-0.312*sin")
for k in [15,25,39,51, 16,26,40,50]:   # odd batch then even batch
    q=Qdiag(k)
    pred=math.log(k*float(b)/(2*math.pi))-amp*math.cos(k*float(b)*float(log2))
    rk=(q-pred)*k
    model=-0.5-float(2*wp/mp.pi)*math.sin(k*float(b)*float(log2))
    par='odd ' if k%2 else 'even'
    print(f" {par} {k:3d}  {q:8.5f}  {rk:+.4f}   {model:+.4f}")

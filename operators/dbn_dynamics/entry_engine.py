#!/usr/bin/env python3
"""§199.3b.2: closed-form entry engine. Q[k,l]=∫_0^{2a} g(t)D_kl(t)dt (+prime part on [log2,2a]).
Uniform elementary integral via COMPLEX exp: ∫ t^p e^{-αt}{cos/sin}(ωt)dt = Re/Im ∫t^p e^{st}dt, s=-α+iω,
∫t^p e^{st}dt = e^{st}·P_p(t,s), P0=1/s, P1=t/s-1/s², P2=t²/s-2t/s²+2/s³.
g-terms (amp·t^{p}·e^{-αt}); D-terms (amp·t^{p}·{cos/sin}(ωt)). Validate 4×4 vs weil_matrices."""
import sys, mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, prime_data
mp.mp.dps=40
a=mp.mpf('0.549'); b=mp.pi/(2*a); log2=mp.log(2); twoa=2*a
Om0=mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi); Phi0=mp.zeta(2,mp.mpf(1)/4); wp=mp.log(2)/mp.sqrt(2)
NLERCH=400
def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3
def II(p,alpha,omega,sc,t0,t1):
    s=mp.mpc(-alpha,omega)
    F=lambda t: mp.e**(s*t)*Pp(p,t,s)
    v=F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
def gterms():
    g=[(-4,0,mp.mpf(-1)/2),(-4,0,mp.mpf(1)/2),(8,0,mp.mpf(0)),(-Om0/2,1,mp.mpf(0)),(-Phi0/4,0,mp.mpf(0))]
    for n in range(NLERCH): g.append((mp.mpf(1)/4/(n+mp.mpf(1)/4)**2,0,2*n+mp.mpf(1)/2))
    return g
def dterms(k,l):
    if k==l:
        return [(2*a*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
    al=k*l*b*(mp.mpf(1)/(k-l)-mp.mpf(1)/(k+l)); ak=-k*l*b*(mp.mpf(1)/(k-l)+mp.mpf(1)/(k+l))
    return [(al,0,l*b,'sin'),(ak,0,k*b,'sin')]
def Q(k,l):
    tot=mp.mpf(0)
    for (ag,pg,alp) in gterms():
        for (ad,pd,om,sc) in dterms(k,l):
            tot+=ag*ad*II(pg+pd,alp,om,sc,mp.mpf(0),twoa)
    # prime part on [log2,2a]: wp*(t-log2) = wp*t - wp*log2
    for (ad,pd,om,sc) in dterms(k,l):
        tot+=wp*ad*II(pd+1,mp.mpf(0),om,sc,log2,twoa) + (-wp*log2)*ad*II(pd,mp.mpf(0),om,sc,log2,twoa)
    return tot
modes=[1,3,5,7]; Qf,Qa,Mass,Vv,xs,w=weil_matrices(a,24,1200)
ev=[c for c in range(24) if (c+1)%2==1][:4]
print(f"§199.3b.2 closed-form entry engine vs mpmath @ a=0.549 (NLERCH={NLERCH})")
maxe=mp.mpf(0)
for ik,k in enumerate(modes):
    for il,l in enumerate(modes):
        if il<ik: continue
        v=Q(k,l); r=Qf[ev[ik],ev[il]]; e=abs(v-r); maxe=max(maxe,e)
        print(f"  Q[{k},{l}] cf={mp.nstr(v,12)} mpmath={mp.nstr(r,12)} err={mp.nstr(e,3)}",flush=True)
print(f"\nmax err={mp.nstr(maxe,3)} (small ⟹ engine CORRECT; residual=Lerch trunc+mpmath grid. Then Arb port.)")

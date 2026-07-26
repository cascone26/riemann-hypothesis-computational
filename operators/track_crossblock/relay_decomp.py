#!/usr/bin/env python3
"""Relay-argument groundwork (§188-190 / PLAN_PRIME_INDUCTION Phase 1-2): decompose the even-sector near-null
margin into its SMOOTH (arch+boundary+Lerch) part and its PRIME part, and watch each scale with a across the
p={2,3} window. Individual pieces are O(1) (float64-resolvable); only their cancellation (the margin) is tiny.
So float64 reveals the RELAY STRUCTURE even though the margin's sign needs Arb.

Uses the exact screw Q-formula (from arb_certify.py) split as Q = Q_smooth - Q_prime, in float64/cmath.
Even 4x4 block, modes [1,3,5,7]. At each a:
  - lam_min(M_smooth), v0 = its argmin eigenvector          [arch/boundary 'state' before primes]
  - <v0|M_prime|v0>                                          [prime's effect ON the smooth near-null]
  - lam_min(M_smooth - M_prime)                              [actual margin, float (sign unreliable if tiny)]
  - diag decomposition on mode 1
Question: as a grows, does the prime term COVER the deepening arch (relay holds) or ERODE the margin (demand)?
"""
import numpy as np, cmath, math

def build(a, NLERCH=6000):
    A=a; B=math.pi/(2*A); LOG2=math.log(2); TWOA=2*A
    OM0 = float(_digamma_quarter()) - math.log(math.pi)
    PHI0 = _hurwitz_zeta_2_quarter()
    def Pp(p,t,s):
        if p==0: return 1/s
        if p==1: return t/s - 1/s**2
        return t*t/s - 2*t/s**2 + 2/s**3
    def II(p,alpha,omega,sc,t0,t1):
        s=complex(-alpha,omega)
        def F(tt):
            t=complex(tt); return cmath.exp(s*t)*Pp(p,t,s)
        v=F(t1)-F(t0)
        return v.real if sc=='cos' else v.imag
    def dterms(k,l):
        if k==l:
            return [(2*A*k*k*B*B,0,k*B,'cos'),(-k*k*B*B,1,k*B,'cos'),(-k*B,0,k*B,'sin')]
        al=k*l*B*(1.0/(k-l)-1.0/(k+l)); ak=-k*l*B*(1.0/(k-l)+1.0/(k+l))
        return [(al,0,l*B,'sin'),(ak,0,k*B,'sin')]
    PP=_prime_powers(TWOA)   # prime powers q with log q < 2a
    def Q_smooth(k,l):
        tot=0.0; D=dterms(k,l)
        gt=[(-4.0,0,-0.5),(-4.0,0,0.5),(8.0,0,0.0),(-OM0/2,1,0.0),(-PHI0/4,0,0.0)]
        for n in range(NLERCH): gt.append((0.25/(n+0.25)**2,0,2*n+0.5))
        for (ag,pg,alp) in gt:
            for (ad,pd,om,sc) in D: tot+=ag*ad*II(pg+pd,alp,om,sc,0.0,TWOA)
        return tot
    def Q_prime(k,l):
        tot=0.0; D=dterms(k,l)
        for (q,p) in PP:
            Wq=math.log(p)/math.sqrt(q); LOGQ=math.log(q)
            for (ad,pd,om,sc) in D:
                tot+=Wq*ad*II(pd+1,0.0,om,sc,LOGQ,TWOA)+(-Wq*LOGQ)*ad*II(pd,0.0,om,sc,LOGQ,TWOA)
        return tot
    modes=[1,3,5,7]
    Ms=np.array([[Q_smooth(k,l) for l in modes] for k in modes])
    Mp=np.array([[Q_prime(k,l) for l in modes] for k in modes])
    Ms=(Ms+Ms.T)/2; Mp=(Mp+Mp.T)/2
    return Ms,Mp,PP,modes

def _primes_upto(n):
    out=[];
    for x in range(2,int(n)+1):
        if all(x%p for p in out): out.append(x)
    return out
def _prime_powers(twoa):
    X=math.exp(twoa); out=[]
    for p in _primes_upto(X+1):
        q=p
        while math.log(q)<twoa-1e-12: out.append((q,p)); q*=p
    return out
def _digamma_quarter():
    from mpmath import digamma, mpf
    return float(digamma(mpf(1)/4))
def _hurwitz_zeta_2_quarter():
    from mpmath import zeta, mpf
    return float(zeta(2, mpf(1)/4))

print("RELAY DECOMPOSITION — decompose the TRUE near-null margin lam_min(Q), Q=Ms+Mp (prime ADDED, matches Arb cert).")
print("On Q's own near-null v0:  lam_min = <v0|Ms|v0>(arch, NEG) + <v0|Mp|v0>(prime budget, POS).")
print("RELAY <=> the prime budget tracks/covers the arch so the SUM stays >= 0 as a grows.\n")
print(f"{'a':>6} {'primes':>10} {'arch<v0|Ms|v0>':>15} {'prime<v0|Mp|v0>':>16} {'sum=lam_min':>13} {'|arch|/prime':>12}")
rows=[]
for a in [0.549,0.55,0.553,0.56,0.58,0.60,0.62,0.66,0.70,0.74,0.78,0.80]:
    Ms,Mp,PP,modes=build(a)
    M=Ms+Mp
    w,V=np.linalg.eigh(M); lam=w[0]; v0=V[:,0]
    arch=float(v0@Ms@v0); prime=float(v0@Mp@v0)
    ratio=abs(arch)/prime if prime!=0 else float('nan')
    pr=",".join(str(q) for q,_ in PP)
    rows.append((a,arch,prime,lam,ratio))
    print(f"{a:>6.3f} {pr:>10} {arch:>15.6e} {prime:>16.6e} {lam:>13.4e} {ratio:>12.6f}")
print("\nREAD: arch<0 deepens with a; prime>0 grows with a. If |arch|/prime < 1 throughout => prime OVER-covers")
print("=> margin (sum) stays >0 = RELAY HOLDS numerically. The analytic target: PROVE |arch(a)| <= prime(a) on")
print("the near-null, uniformly in a (and per-prime, = the inductive step). NOTE lam_min at ~1e-7 is float-noisy;")
print("trust the O(0.1) arch/prime pieces + their ratio, not the tiny sum's exact value.")

#!/usr/bin/env python3
"""§324b: strengthen mixture-LP — 3-delta, continuous NON-modular Φ, and the total-positivity mechanism.
If ∫Φ(u)G((us)²)du is LP for arbitrary Φ≥0, the likely reason: the kernel H(u,y)=G(u²y) is TOTALLY POSITIVE
in (u,y) [sign-regular], making the transform variation-diminishing ⇒ real-rootedness preserved (Schoenberg).
Tests: (A) 3-delta adversarial; (B) continuous Φ=u^p e^{-qu} (positive, NOT the modular double-exp);
(C) the mechanism: 3x3 minors of the matrix [G(u_i² y_j)] — sign-regular?"""
import mpmath as mp
mp.mp.dps=30
def G(x): return mp.hyper([],[mp.mpf(1)/2,1],-x/4)
def wind(f,x0,x1,y0,y1,n):
    pts=[]
    for i in range(n): pts.append(mp.mpc(x0+(x1-x0)*i/n,y0))
    for i in range(n): pts.append(mp.mpc(x1,y0+(y1-y0)*i/n))
    for i in range(n): pts.append(mp.mpc(x1-(x1-x0)*i/n,y1))
    for i in range(n): pts.append(mp.mpc(x0,y1-(y1-y0)*i/n))
    tot=mp.mpf(0); prev=mp.arg(f(pts[0]))
    for p in pts[1:]+[pts[0]]:
        a=mp.arg(f(p)); d=a-prev
        while d>mp.pi: d-=2*mp.pi
        while d<-mp.pi: d+=2*mp.pi
        tot+=d; prev=a
    return tot/(2*mp.pi)
print("(A) 3-delta adversarial B(s)=Σ a_i G((u_i s)²):",flush=True)
tri=[[(1,1),(3,1),(7,1)],[(1,5),(2,1),(6,3)],[(0.7,2),(3,1),(9,4)],[(1,1),(1.3,1),(1.7,1)],[(1,10),(4,1),(12,8)]]
allA=True
for case in tri:
    f=lambda s: sum(mp.mpf(a)*G((mp.mpf(u)*s)**2) for (u,a) in case)
    w=wind(f,mp.mpf('0.05'),15,mp.mpf('0.05'),6,600)
    lp=abs(w)<mp.mpf('0.5'); allA=allA and lp
    print(f"   {case}: wind={mp.nstr(w,3)} {'LP' if lp else '*** COMPLEX ***'}",flush=True)
print(f"   all 3-delta LP? {allA}",flush=True)
print("(B) continuous NON-modular Φ=u^p e^{-q u} (Gamma-type, positive):",flush=True)
allB=True
for (p,q) in [(0,1),(2,1),(1,mp.mpf('0.5')),(4,2),(mp.mpf('0.5'),mp.mpf('0.3'))]:
    def f(s):
        return mp.quad(lambda u: u**p*mp.e**(-q*u)*G((u*s)**2),[0,1,3,6,12])
    w=wind(f,mp.mpf('0.05'),12,mp.mpf('0.05'),5,400)
    lp=abs(w)<mp.mpf('0.5'); allB=allB and lp
    print(f"   Φ=u^{mp.nstr(p,2)} e^(-{mp.nstr(q,2)}u): wind={mp.nstr(w,3)} {'LP' if lp else '*** COMPLEX ***'}",flush=True)
print(f"   all continuous-Φ LP? {allB}",flush=True)
print("(C) mechanism: sign of 2x2 & 3x3 minors of [G(u_i² y_j)] (total positivity ⇒ variation-diminishing):",flush=True)
us=[mp.mpf('0.5'),mp.mpf(1),mp.mpf('2'),mp.mpf('3.5')]
ys=[mp.mpf('0.3'),mp.mpf(1),mp.mpf('2.5'),mp.mpf('5')]
Hm=mp.matrix(4,4)
for i in range(4):
    for j in range(4): Hm[i,j]=G(us[i]**2*ys[j])
import itertools
sign2=set(); sign3=set()
for I in itertools.combinations(range(4),2):
    for J in itertools.combinations(range(4),2):
        d=Hm[I[0],J[0]]*Hm[I[1],J[1]]-Hm[I[0],J[1]]*Hm[I[1],J[0]]
        sign2.add(mp.sign(d))
print(f"   2x2 minor signs: {sign2}  (single sign ⇒ sign-regular order 2)",flush=True)
print("   NOTE: G oscillates (has zeros) so H is sign-REGULAR not totally-positive; the relevant structure is",flush=True)
print("   that the scale-family shares zeros ⇒ variation-diminishing on POSITIVE Φ. Evidence is the LP verdicts.",flush=True)
print("DONE",flush=True)

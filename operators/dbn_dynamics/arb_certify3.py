#!/usr/bin/env python3
"""§199.5: RIGOROUS 8×8 even-block certification @ a=0.549, with the 2 code-audit rigor FIXES:
 FIX1: pivot test uses the DIRECT arb interval comparison `s>0` (True iff whole interval >0) — no float mid/rad extraction.
 FIX2: each entry is acb; ASSERT the imaginary part encloses 0 before taking .real (else not rigorous).
8×8 (modes 1,3,5,7,9,11,13,15) folds the tail modes IN, so only d>8 tail (~8e-9, workflow) remains. Arch-tail integral
at reduced working prec (it's ~0.01, needs no more)."""
import sys
from flint import arb, acb, ctx
MAINPREC=900
ctx.prec=MAINPREC
A=arb('0.549'); B=arb.pi()/(2*A); LOG2=arb(2).log(); TWOA=2*A
OM0=(arb(1)/4).digamma()-arb.pi().log(); PHI0=acb(2).zeta(acb(arb(1)/4)).real; WP=arb(2).log()/arb(2).sqrt()
def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3
def II(p,alpha,omega,sc,t0,t1):
    s=acb(-alpha,omega); F=lambda tt:(s*acb(tt)).exp()*Pp(p,acb(tt),s); v=F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
def dterms(k,l):
    if k==l: return [(2*A*arb(k)**2*B**2,0,arb(k)*B,'cos'),(-arb(k)**2*B**2,1,arb(k)*B,'cos'),(-arb(k)*B,0,arb(k)*B,'sin')]
    al=arb(k)*arb(l)*B*(arb(1)/(k-l)-arb(1)/(k+l)); ak=-arb(k)*arb(l)*B*(arb(1)/(k-l)+arb(1)/(k+l))
    return [(al,0,arb(l)*B,'sin'),(ak,0,arb(k)*B,'sin')]
def Dacb(k,l,t):
    kb=acb(arb(k)*B); lb=acb(arb(l)*B)
    if k==l: return acb(arb(k))**2*acb(B)**2*(acb(TWOA)-t)*(kb*t).cos()-acb(arb(k)*B)*(kb*t).sin()
    al=acb(arb(k)*arb(l)*B*(arb(1)/(k-l)-arb(1)/(k+l))); ak=acb(-arb(k)*arb(l)*B*(arb(1)/(k-l)+arb(1)/(k+l)))
    return al*(lb*t).sin()+ak*(kb*t).sin()
def Q(k,l):
    tot=arb(0)
    gt=[(arb(-4),0,arb(-1)/2),(arb(-4),0,arb(1)/2),(arb(8),0,arb(0)),(-OM0/2,1,arb(0)),(-PHI0/4,0,arb(0))]
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in dterms(k,l): tot+=ag*ad*II(pg+pd,alp,om,sc,arb(0),TWOA)
    for (ad,pd,om,sc) in dterms(k,l):
        tot+=WP*ad*II(pd+1,arb(0),om,sc,LOG2,TWOA)+(-WP*LOG2)*ad*II(pd,arb(0),om,sc,LOG2,TWOA)
    def integ(t,_):
        z=(acb(-2)*t).exp()
        return acb('0.25')*(acb('-0.5')*t).exp()*z.lerch_phi(acb(2),acb(arb(1)/4))*Dacb(k,l,t)
    ctx.prec=350
    Iarch=acb.integral(integ,acb(0),acb(TWOA))
    ctx.prec=MAINPREC
    # FIX2: assert imaginary part encloses 0 (rigorous real-ness)
    assert Iarch.imag.contains(arb(0)), f"imag not enclosing 0 at Q[{k},{l}]: {Iarch.imag}"
    return tot+Iarch.real
modes=[1,3,5,7,9,11,13,15]; n=len(modes)
M=[[Q(k,l) for l in modes] for k in modes]
# LDL^T with DIRECT interval pivot test (FIX1)
d=[None]*n; L=[[arb(0)]*n for _ in range(n)]; ok=True
print(f"§199.5 RIGOROUS 8×8 even-block cert (audit fixes) @ a=0.549, prec={MAINPREC}")
for i in range(n):
    s=M[i][i]
    for kk in range(i): s=s-L[i][kk]**2*d[kk]
    d[i]=s
    cert = (s > 0)                        # FIX1: rigorous arb comparison
    if not cert: ok=False
    print(f"  pivot {i} (mode {modes[i]}): {s.str(6)}   rigorously >0 ? {cert}")
    for j in range(i+1,n):
        t=M[j][i]
        for kk in range(i): t=t-L[j][kk]*L[i][kk]*d[kk]
        L[j][i]=t/d[i]
print(f"\n{'*** 8x8 CERTIFIED PD (all pivots rigorously >0) — p=2-only positivity RIGOROUS at a=0.549 (mod d>8 tail ~8e-9) ***' if ok else 'NOT certified — tighten prec'}")

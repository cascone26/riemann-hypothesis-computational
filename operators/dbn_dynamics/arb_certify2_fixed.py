#!/usr/bin/env python3
"""§199.3c v2: RIGOROUS Arb certification, even 4×4 @ a=0.549, using NATIVE acb.lerch_phi + acb.integral for the
smooth arch-tail (kills the loose Lerch-tail-bound that made pivots too wide). Elementary parts closed-form (II)."""
import sys
from flint import arb, acb, ctx
ctx.prec = 500
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
def Dacb(k,l,t):   # D_kl(t) as acb (t real acb)
    kb=acb(arb(k)*B); lb=acb(arb(l)*B)
    if k==l:
        return acb(arb(k))**2*acb(B)**2*(acb(TWOA)-t)*(kb*t).cos() - acb(arb(k)*B)*(kb*t).sin()
    al=acb(arb(k)*arb(l)*B*(arb(1)/(k-l)-arb(1)/(k+l))); ak=acb(-arb(k)*arb(l)*B*(arb(1)/(k-l)+arb(1)/(k+l)))
    return al*(lb*t).sin()+ak*(kb*t).sin()
def Q(k,l):
    tot=arb(0)
    # elementary g: boundary + arch_core + (-Phi0/4 const) + prime
    gt=[(arb(-4),0,arb(-1)/2),(arb(-4),0,arb(1)/2),(arb(8),0,arb(0)),(-OM0/2,1,arb(0)),(-PHI0/4,0,arb(0))]
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in dterms(k,l): tot+=ag*ad*II(pg+pd,alp,om,sc,arb(0),TWOA)
    for (ad,pd,om,sc) in dterms(k,l):
        tot+=WP*ad*II(pd+1,arb(0),om,sc,LOG2,TWOA)+(-WP*LOG2)*ad*II(pd,arb(0),om,sc,LOG2,TWOA)
    # arch-tail Lerch: ¼∫_0^{2a} e^{-t/2} Φ(e^{-2t},2,¼) D(t) dt  (rigorous, native lerch_phi)
    def integ(t,_):
        z=(acb(-2)*t).exp()
        return acb('0.25')*(acb('-0.5')*t).exp()*z.lerch_phi(acb(2),acb(arb(1)/4))*Dacb(k,l,t)
    ctx.prec=350
    Ia=acb.integral(integ,acb(0),acb(TWOA)); ctx.prec=500
    assert Ia.imag.contains(arb(0)), 'imag not 0'
    tot+=Ia.real
    return tot
modes=[1,3,5,7]; M=[[Q(k,l) for l in modes] for k in modes]
d=[None]*4; L=[[arb(0)]*4 for _ in range(4)]; ok=True; piv=[]
for i in range(4):
    s=M[i][i]
    for kk in range(i): s=s-L[i][kk]**2*d[kk]
    d[i]=s; cert=(s>0); piv.append((s.str(6),cert))
    if not cert: ok=False
    for j in range(i+1,4):
        t=M[j][i]
        for kk in range(i): t=t-L[j][kk]*L[i][kk]*d[kk]
        L[j][i]=t/d[i]
print("§199.3c v2 RIGOROUS Arb cert (native lerch_phi + acb.integral), even 4×4 @ a=0.549, prec=500")
for i,(sv,c) in enumerate(piv): print(f"  pivot {i}: {sv}  rigorously>0? {c}")
print(f"\n{'*** CERTIFIED PD — p=2-only positivity RIGOROUS at a=0.549 ***' if ok else 'not certified (tighten prec/integral tol)'}")

#!/usr/bin/env python3
"""§199.3c v2-opt2: MOST AGGRESSIVE optimization - split precision contexts.
Elementary parts at 500 bits (II calculations).
Lerch integral at 250 bits global context (sufficient for ~0.01 contribution).
This should cut integration time by ~70-80%."""
import sys
from flint import arb, acb, ctx
import time

start_time = time.time()

# Global high precision for closed-form parts
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

def Dacb_high(k,l,t):   # D_kl(t) as acb at HIGH precision
    kb=acb(arb(k)*B); lb=acb(arb(l)*B)
    if k==l:
        return acb(arb(k))**2*acb(B)**2*(acb(TWOA)-t)*(kb*t).cos() - acb(arb(k)*B)*(kb*t).sin()
    al=acb(arb(k)*arb(l)*B*(arb(1)/(k-l)-arb(1)/(k+l))); ak=acb(-arb(k)*arb(l)*B*(arb(1)/(k-l)+arb(1)/(k+l)))
    return al*(lb*t).sin()+ak*(kb*t).sin()

def Q(k,l):
    tot=arb(0)
    # elementary g: boundary + arch_core + (-Phi0/4 const) + prime
    # All at high precision (500 bits)
    gt=[(arb(-4),0,arb(-1)/2),(arb(-4),0,arb(1)/2),(arb(8),0,arb(0)),(-OM0/2,1,arb(0)),(-PHI0/4,0,arb(0))]
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in dterms(k,l): tot+=ag*ad*II(pg+pd,alp,om,sc,arb(0),TWOA)
    for (ad,pd,om,sc) in dterms(k,l):
        tot+=WP*ad*II(pd+1,arb(0),om,sc,LOG2,TWOA)+(-WP*LOG2)*ad*II(pd,arb(0),om,sc,LOG2,TWOA)

    # arch-tail Lerch: ¼∫_0^{2a} e^{-t/2} Φ(e^{-2t},2,¼) D(t) dt
    # AGGRESSIVE OPTIMIZATION: switch to 250-bit context for the integral
    # This is safe because the contribution is small (~0.01) and the integrand is smooth
    saved_prec = ctx.prec
    ctx.prec = 250

    # Pre-compute constants at low precision (they're reused in integ calls)
    A_low = arb('0.549'); B_low = arb.pi()/(2*A_low); TWOA_low = 2*A_low

    # Recompute D at low precision for use in integral
    def integ(t,_):
        kb_low=acb(arb(k)*B_low); lb_low=acb(arb(l)*B_low)
        if k==l:
            D = acb(arb(k))**2*acb(B_low)**2*(acb(TWOA_low)-t)*(kb_low*t).cos() - acb(arb(k)*B_low)*(kb_low*t).sin()
        else:
            al_low=acb(arb(k)*arb(l)*B_low*(arb(1)/(k-l)-arb(1)/(k+l)))
            ak_low=acb(-arb(k)*arb(l)*B_low*(arb(1)/(k-l)+arb(1)/(k+l)))
            D = al_low*(lb_low*t).sin()+ak_low*(kb_low*t).sin()
        z=(acb(-2)*t).exp()
        return acb('0.25')*(acb('-0.5')*t).exp()*z.lerch_phi(acb(2),acb(arb(1)/4))*D

    lerch_contrib = acb.integral(integ, acb(0), acb(TWOA_low)).real
    ctx.prec = saved_prec  # Restore to 500

    tot += arb(str(lerch_contrib))
    return tot

modes=[1,3,5,7]; M=[[Q(k,l) for l in modes] for k in modes]
d=[None]*4; L=[[arb(0)]*4 for _ in range(4)]; ok=True; piv=[]

for i in range(4):
    s=M[i][i]
    for kk in range(i): s=s-L[i][kk]**2*d[kk]
    d[i]=s; lo=arb(s.mid())-arb(s.rad()); piv.append((float(s.mid()),float(s.rad())))
    if not(lo>0): ok=False
    for j in range(i+1,4):
        t=M[j][i]
        for kk in range(i): t=t-L[j][kk]*L[i][kk]*d[kk]
        L[j][i]=t/d[i]

elapsed = time.time() - start_time

print("§199.3c v2-opt2 AGGRESSIVE-OPT Arb cert (prec-switched integral), even 4×4 @ a=0.549")
print(f"Elementary prec=500, integral prec=250, runtime={elapsed:.1f}s")
for i,(m,r) in enumerate(piv):
    status = '>0 OK' if m-r>0 else 'FAIL'
    print(f"  pivot {i}: {m:.6e} +/- {r:.2e}  lower={m-r:.6e}  {status}")
print(f"\n{'*** CERTIFIED PD — p=2-only positivity RIGOROUS at a=0.549 ***' if ok else 'not certified (tighten prec/integral tol)'}")

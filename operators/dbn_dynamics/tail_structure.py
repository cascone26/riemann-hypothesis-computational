#!/usr/bin/env python3
"""§201: empirical tail structure for the tail lemma. Compute even-sector Q[k,k] (diagonal) and off-diagonal row
sum Σ_{l≠k}|Q[k,l]| for high odd k, via the fast closed-form entry engine. Does the diagonal DOMINATE (Gershgorin
disc excludes 0) for k>K0? That + growth (arch symbol Ω(r)~log r) = the tail lemma. a=0.549 (window endpoint, tightest)."""
import sys, mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
mp.mp.dps=30
a=mp.mpf('0.549'); b=mp.pi/(2*a); log2=mp.log(2); twoa=2*a
Om0=mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi); Phi0=mp.zeta(2,mp.mpf(1)/4); wp=mp.log(2)/mp.sqrt(2); NL=3000
def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3
def II(p,al,om,sc,t0,t1):
    s=mp.mpc(-al,om); F=lambda tt: mp.e**(s*tt)*Pp(p,tt,s); v=F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
def dterms(k,l):
    if k==l: return [(2*a*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
    al=k*l*b*(mp.mpf(1)/(k-l)-mp.mpf(1)/(k+l)); ak=-k*l*b*(mp.mpf(1)/(k-l)+mp.mpf(1)/(k+l))
    return [(al,0,l*b,'sin'),(ak,0,k*b,'sin')]
def Q(k,l):
    tot=mp.mpf(0)
    gt=[(-4,0,mp.mpf(-1)/2),(-4,0,mp.mpf(1)/2),(8,0,mp.mpf(0)),(-Om0/2,1,mp.mpf(0)),(-Phi0/4,0,mp.mpf(0))]
    for n in range(NL): gt.append((mp.mpf(1)/4/(n+mp.mpf(1)/4)**2,0,2*n+mp.mpf(1)/2))
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in dterms(k,l): tot+=ag*ad*II(pg+pd,alp,om,sc,mp.mpf(0),twoa)
    for (ad,pd,om,sc) in dterms(k,l):
        tot+=wp*ad*II(pd+1,mp.mpf(0),om,sc,log2,twoa)+(-wp*log2)*ad*II(pd,mp.mpf(0),om,sc,log2,twoa)
    return tot
odds=list(range(1,60,2))
print(f"§201 tail structure @ a=0.549 (odd modes)   Mass=a*I so whiten: Q[k,k]/a etc.")
print(f"{'k':>4} {'Q[k,k]':>13} {'Q[k,k]/k^2':>12} {'offdiag rowsum':>15} {'Gershgorin lo':>14} {'diag/rowsum':>11}")
for k in [1,3,5,7,11,15,21,31,41,51]:
    qkk=Q(k,k); rs=mp.fsum(abs(Q(k,l)) for l in odds if l!=k)
    print(f"{k:>4} {mp.nstr(qkk,6):>13} {mp.nstr(qkk/k**2,5):>12} {mp.nstr(rs,6):>15} {mp.nstr(qkk-rs,4):>14} {mp.nstr(qkk/rs,4):>11}",flush=True)
print("\nREAD: if Q[k,k]-rowsum >0 (Gershgorin lo >0) for k>=K0, the tail block modes>K0 are PD by Gershgorin ⟹ tail")
print("lemma. Q[k,k]/k^2 -> const? or ~log k? = the growth law. diag/rowsum>1 = strict diagonal dominance.")

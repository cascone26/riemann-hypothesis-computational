#!/usr/bin/env python3
"""piece-5 industrialization: generalize the screw Weil cert (operators/dbn_dynamics/arb_certify.py,
which is arch+Lerch+p=2 ONLY) to ALL prime powers q with log q < 2a. New file — canonical cert untouched.

The p=2 prime term in arb_certify.py is:
    WP = log2/sqrt(2)             # = Lambda(2)/sqrt(2)
    for (ad,pd,om,sc) in D:
        tot += WP*ad*II(pd+1,0,om,sc,LOG2,TWOA) + (-WP*LOG2)*ad*II(pd,0,om,sc,LOG2,TWOA)
i.e. weight Lambda(2)/sqrt(2), position log2, screw ramp (t-log2) over [log2, 2a]. The structural
generalization: sum the SAME term over every prime power q=p^m with log q < 2a, weight Lambda(q)/sqrt(q)
= log(p)/sqrt(q), position log q, integrated over [log q, 2a].

CORRECTNESS GATE: at a=0.549 (2a=log3), only q=2 has log q<2a (q=3 sits exactly at 2a, zero-width) — so this
MUST reproduce arb_certify.py's p=2-only pivots. Verified before trusting any a>1/2 log 3 extension.

Usage: arb_certify_multiprime.py [--a 0.549] [--parity even|odd]
"""
import argparse
from flint import arb, acb, ctx
ctx.prec = 1400

p = argparse.ArgumentParser()
p.add_argument('--a', type=str, default='0.549')
p.add_argument('--parity', choices=['even','odd'], default='even')
p.add_argument('--prec', type=int, default=1400)
p.add_argument('--nlerch', type=int, default=8000)
p.add_argument('--k0', type=int, default=0)
args = p.parse_args()
ctx.prec = args.prec

A = arb(args.a); B = arb.pi()/(2*A); LOG2 = arb(2).log(); TWOA = 2*A
OM0 = (arb(1)/4).digamma() - arb.pi().log()
PHI0 = acb(2).zeta(acb(arb(1)/4)).real
NLERCH = args.nlerch

def prime_powers_upto_log(twoa_float):
    """(q, logp) for prime powers q=p^m with log(q) < 2a (strict)."""
    import math
    X = math.exp(twoa_float)
    primes = []; x = 2
    while x <= X + 1:
        if all(x % pr for pr in primes): primes.append(x)
        x += 1
    out = []
    for pr in primes:
        q = pr
        while math.log(q) < twoa_float - 1e-12:   # strict: exclude q with log q == 2a (zero-width)
            out.append((q, pr)); q *= pr
    return out

PP = prime_powers_upto_log(2*float(args.a))
# weights Lambda(q)/sqrt(q) = log(p)/sqrt(q); positions log q
PRIMES = [(arb(pr).log()/arb(q).sqrt(), arb(q).log()) for (q, pr) in PP]

def Pp(p_,t,s):
    if p_==0: return 1/s
    if p_==1: return t/s - 1/s**2
    return t*t/s - 2*t/s**2 + 2/s**3
def II(p_, alpha, omega, sc, t0, t1):
    s = acb(-alpha, omega)
    def F(tt):
        t = acb(tt); return (s*t).exp()*Pp(p_,t,s)
    v = F(t1) - F(t0)
    return v.real if sc=='cos' else v.imag
def dterms(k,l):
    kb=arb(k)*B; lb=arb(l)*B
    if k==l:
        return [(2*A*arb(k)**2*B**2,0,kb,'cos'),(-arb(k)**2*B**2,1,kb,'cos'),(-arb(k)*B,0,kb,'sin')]
    al=arb(k)*arb(l)*B*(arb(1)/(k-l)-arb(1)/(k+l)); ak=-arb(k)*arb(l)*B*(arb(1)/(k-l)+arb(1)/(k+l))
    return [(al,0,arb(l)*B,'sin'),(ak,0,arb(k)*B,'sin')]
def Dmax(k,l):
    tot=arb(0)
    for (ad,pd,om,sc) in dterms(k,l): tot += abs(ad)*(TWOA if pd==1 else arb(1))
    return tot
def Q(k,l):
    tot=arb(0); D=dterms(k,l)
    gt=[(arb(-4),0,arb(-1)/2),(arb(-4),0,arb(1)/2),(arb(8),0,arb(0)),(-OM0/2,1,arb(0)),(-PHI0/4,0,arb(0))]
    for n in range(NLERCH): gt.append((arb(1)/4/(arb(n)+arb(1)/4)**2,0,2*arb(n)+arb(1)/2))
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in D: tot += ag*ad*II(pg+pd,alp,om,sc,arb(0),TWOA)
    N=NLERCH; dm=Dmax(k,l)
    tailbound = dm/8 * (arb(1)/(2*(arb(N)-1)**2))
    tot += tailbound*arb(0,1)
    # GENERALIZED prime part: sum over all prime powers q with log q < 2a
    for (Wq, LOGQ) in PRIMES:
        for (ad,pd,om,sc) in D:
            tot += Wq*ad*II(pd+1,arb(0),om,sc,LOGQ,TWOA) + (-Wq*LOGQ)*ad*II(pd,arb(0),om,sc,LOGQ,TWOA)
    return tot

if args.parity=='even':
    modes = list(range(2*args.k0-1, 0, -2)) if args.k0>0 else [5,7,3,1]  # odd modes 1..2k0-1, large-first
else:
    modes = list(range(2*args.k0-1, 0, -2)) if args.k0>0 else [1]
n_modes=len(modes)
M=[[Q(k,l) for l in modes] for k in modes]
d=[None]*n_modes; L=[[arb(0)]*n_modes for _ in range(n_modes)]; ok=True; piv=[]
for i in range(n_modes):
    s=M[i][i]
    for kk in range(i): s=s-L[i][kk]**2*d[kk]
    d[i]=s; c=(s>0); piv.append((s.str(6),c))
    if not c: ok=False
    for j in range(i+1,n_modes):
        t=M[j][i]
        for kk in range(i): t=t-L[j][kk]*L[i][kk]*d[kk]
        L[j][i]=t/d[i]
primestr=", ".join(f"{q}" for (q,pr) in PP) or "(none)"
print(f"MULTIPRIME screw cert, {args.parity.upper()} sector, a={args.a} (2a={2*float(args.a):.5f}), prec=1400, NLERCH={NLERCH}")
print(f"prime powers included (log q < 2a): {primestr}   [{len(PP)} terms]")
for i,(sv,c) in enumerate(piv):
    print(f"  pivot {i} (mode {modes[i]}): {sv}  rigorously>0? {c}")
print(f"\n{'*** CERTIFIED PD ***' if ok else 'NOT certified'}")

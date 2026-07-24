#!/usr/bin/env python3
"""§199.3c/§200: RIGOROUS Arb certification of Weil blocks PD at a=0.549 (even or odd parity).
Entries via the §199.3b.2 engine ported to acb intervals + rigorous Lerch tail bound. Mass=a·I (orthogonal
sines). Certifies Q_even (even parity, modes 5,7,3,1) or Q_odd (odd parity, mode 1) to be PD, proving
p=2-only positivity at a=0.549 for the corresponding sector.
Usage: arb_certify.py [--parity even|odd]  (default: even)"""
import sys
import argparse
from flint import arb, acb, ctx
ctx.prec = 1400

# Parse parity argument
parser = argparse.ArgumentParser(description='Rigorous Arb certification of Weil blocks')
parser.add_argument('--parity', choices=['even', 'odd'], default='even', help='Parity sector to certify')
args = parser.parse_args()
parity = args.parity

A = arb('0.549'); B = arb.pi()/(2*A); LOG2 = arb(2).log(); TWOA = 2*A
OM0 = (arb(1)/4).digamma() - arb.pi().log()
PHI0 = acb(2).zeta(acb(arb(1)/4)).real
WP = arb(2).log()/arb(2).sqrt()
NLERCH = 8000
def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s - 1/s**2
    return t*t/s - 2*t/s**2 + 2/s**3
def II(p, alpha, omega, sc, t0, t1):   # ∫_{t0}^{t1} t^p e^{-alpha t}{cos/sin}(omega t) dt  (arb)
    s = acb(-alpha, omega)
    def F(tt):
        t = acb(tt)
        return (s*t).exp()*Pp(p,t,s)
    v = F(t1) - F(t0)
    return v.real if sc=='cos' else v.imag
def dterms(k,l):
    kb=arb(k)*B; lb=arb(l)*B
    if k==l:
        return [(2*A*arb(k)**2*B**2,0,kb,'cos'),(-arb(k)**2*B**2,1,kb,'cos'),(-arb(k)*B,0,kb,'sin')]
    al=arb(k)*arb(l)*B*(arb(1)/(k-l)-arb(1)/(k+l)); ak=-arb(k)*arb(l)*B*(arb(1)/(k-l)+arb(1)/(k+l))
    return [(al,0,arb(l)*B,'sin'),(ak,0,arb(k)*B,'sin')]
def Dmax(k,l):   # bound max|D_kl(t)| on [0,2a]: Σ|amp|·(2a if p=1 else 1)
    tot=arb(0)
    for (ad,pd,om,sc) in dterms(k,l): tot += abs(ad)*(TWOA if pd==1 else arb(1))
    return tot
def Q(k,l):
    tot=arb(0); D=dterms(k,l)
    # g smooth part: boundary + arch_core + (-Phi0/4) + Lerch(NLERCH terms)
    gt=[(arb(-4),0,arb(-1)/2),(arb(-4),0,arb(1)/2),(arb(8),0,arb(0)),(-OM0/2,1,arb(0)),(-PHI0/4,0,arb(0))]
    for n in range(NLERCH): gt.append((arb(1)/4/(arb(n)+arb(1)/4)**2,0,2*arb(n)+arb(1)/2))
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in D: tot += ag*ad*II(pg+pd,alp,om,sc,arb(0),TWOA)
    # Lerch rigorous tail (n>=NLERCH): |term| <= (1/4)(n+1/4)^-2 * Dmax/(2n+1/2); Σ_{n>=N} <= (Dmax/8)Σ n^-3
    N=NLERCH; dm=Dmax(k,l)
    tailbound = dm/8 * (arb(1)/(2*(arb(N)-1)**2))     # Σ_{n>=N} n^-3 <= 1/(2(N-1)^2)
    tot += tailbound*arb(0,1)   # [-tailbound, tailbound]: arb(mid=0,rad=1)=[-1,1]  # enclose ±tailbound
    # prime part on [log2,2a]
    for (ad,pd,om,sc) in D:
        tot += WP*ad*II(pd+1,arb(0),om,sc,LOG2,TWOA) + (-WP*LOG2)*ad*II(pd,arb(0),om,sc,LOG2,TWOA)
    return tot
# Select modes by parity (large-diagonal first to avoid LDL^T interval amplification)
if parity == 'even':
    modes = [5, 7, 3, 1]  # even parity sector (odd sine indices, large diag first)
else:
    modes = [1]  # odd parity sector (single negative arch mode, generous margin per LEDGER §200)

n_modes = len(modes)
M = [[Q(k, l) for l in modes] for k in modes]

# LDL^T pivots (arb), certify lower bound >0
d = [None] * n_modes
L = [[arb(0)] * n_modes for _ in range(n_modes)]
ok = True
piv = []
for i in range(n_modes):
    s = M[i][i]
    for kk in range(i):
        s = s - L[i][kk]**2 * d[kk]
    d[i] = s
    cert = (s > 0)
    piv.append((s.str(6), cert))
    if not cert:
        ok = False
    for j in range(i + 1, n_modes):
        t = M[j][i]
        for kk in range(i):
            t = t - L[j][kk] * L[i][kk] * d[kk]
        L[j][i] = t / d[i]

print(f"§199.3c/§200 RIGOROUS Arb certification, {parity.upper()} sector @ a=0.549, prec=1400, NLERCH=8000, modes={modes}")
for i, (sv, c) in enumerate(piv):
    print(f"  pivot {i} (mode {modes[i]}): {sv}  rigorously>0? {c}")
print(f"\n{'*** CERTIFIED PD (all pivots lower-bound >0) — p=2-only positivity RIGOROUS at a=0.549 ***' if ok else 'not certified'}")

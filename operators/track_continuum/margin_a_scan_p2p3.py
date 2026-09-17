#!/usr/bin/env python3
"""
cycle 343 — p={2,3} even-sector margin(a) fine scan for MONOTONICITY (the interval-theorem route, cy342).

Copies sos_affine_compose._mk_block (p=2-only) and ADDS the p=3 screw term (weight log3/sqrt3, ramp [log3,2a]).
Interval (0.549,0.567) is SMOOTH (next threshold ½log4=0.693), so margin(a) is a smooth curve. Goal: confirm
d(margin)/da < 0 (monotone decreasing) with NO local bump — then the Arb-certified endpoint a=0.567 (cy341)
would cover the whole interval (odd sector comfortable). REPORT-tier float (mpmath dps30) for the SIGN/shape;
rigorous d(margin)/da<=0 bound is the follow-up if monotone.
"""
import mpmath as mp
mp.mp.dps = 30

def mk_block_p2p3(a, NL=4000):
    a = mp.mpf(a); b = mp.pi/(2*a); log2 = mp.log(2); log3 = mp.log(3); twoa = 2*a
    w2 = log2/mp.sqrt(2); w3 = log3/mp.sqrt(3)
    Om0 = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi); Phi0 = mp.zeta(2, mp.mpf(1)/4)
    def Pp(p,t,s):
        if p==0: return 1/s
        if p==1: return t/s - 1/s**2
        return t*t/s - 2*t/s**2 + 2/s**3
    def II(p,al,om,sc,t0,t1):
        s = mp.mpc(-al,om); F = lambda tt: mp.e**(s*tt)*Pp(p,tt,s); v = F(t1)-F(t0)
        return v.real if sc=='cos' else v.imag
    def dterms(k,l):
        if k==l: return [(2*a*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
        al = k*l*b*(mp.mpf(1)/(k-l)-mp.mpf(1)/(k+l)); ak = -k*l*b*(mp.mpf(1)/(k-l)+mp.mpf(1)/(k+l))
        return [(al,0,l*b,'sin'),(ak,0,k*b,'sin')]
    def Q(k,l):
        tot = mp.mpf(0); D = dterms(k,l)
        gt = [(-4,0,mp.mpf(-1)/2),(-4,0,mp.mpf(1)/2),(8,0,mp.mpf(0)),
              (-Om0/2,1,mp.mpf(0)),(-Phi0/4,0,mp.mpf(0))]
        for n in range(NL): gt.append((mp.mpf(1)/4/(n+mp.mpf(1)/4)**2,0,2*n+mp.mpf(1)/2))
        for (ag,pg,alp) in gt:
            for (ad,pd,om,sc) in D: tot += ag*ad*II(pg+pd,alp,om,sc,0,twoa)
        # p=2 screw ramp [log2,2a]
        for (ad,pd,om,sc) in D:
            tot += w2*ad*II(pd+1,0,om,sc,log2,twoa) + (-w2*log2)*ad*II(pd,0,om,sc,log2,twoa)
        # p=3 screw ramp [log3,2a]  (active when 2a>log3, i.e. a>0.5493)
        if twoa > log3:
            for (ad,pd,om,sc) in D:
                tot += w3*ad*II(pd+1,0,om,sc,log3,twoa) + (-w3*log3)*ad*II(pd,0,om,sc,log3,twoa)
        return tot
    return Q

modes = [1,3,5,7,9,11]   # even sector (odd sine indices), K0=6
def margin(a):
    Q = mk_block_p2p3(a)
    M = mp.matrix(len(modes),len(modes))
    for i,k in enumerate(modes):
        for j,l in enumerate(modes):
            M[i,j] = Q(k,l)
    M = (M + M.T)/2
    ev = mp.eigsy(M, eigvals_only=True)
    return min(ev)

print("p={2,3} even-sector margin(a) = smallest eigenvalue (modes 1..11), mpmath dps30:")
print(f"{'a':>8} {'margin (smallest eig)':>24}")
avals = [mp.mpf(x)/10000 for x in range(5490, 5675, 30)]  # 0.5490..0.5670 step 0.003
prev = None; mono = True
for a in avals:
    m = margin(a)
    flag = ''
    if prev is not None and m > prev: flag = '  <-- INCREASE (breaks monotonicity!)'; mono = False
    print(f"{float(a):>8.4f} {mp.nstr(m,8):>24}{flag}")
    prev = m
print(f"\nmonotone decreasing on [0.549,0.567]? {mono}")
print("If YES + a=0.567 Arb endpoint cert (cy341) ⟹ even-sector positivity on the FULL interval (interval theorem).")

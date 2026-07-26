#!/usr/bin/env python3
"""§324: the DECISIVE mixture-LP experiment. B_Φ(s)=∫Φ(u)G((su)²)du, G=₀F₂(;½,1;−x/4)∈LP.
Is the mixture LP for EVERY Φ≥0 (⇒ unconditional, RH-free) or only special Φ (⇒ RH-adjacent)?
Test adversarial positive Φ: two-point masses at (u1,u2) with varied weights/separation — the extreme
case of a positive measure. If two-delta mixtures a·G((u1 s)²)+b·G((u2 s)²) can have COMPLEX zeros,
the mixture-LP lemma is FALSE ⇒ B's real-rootedness uses the modular Φ (not free). Check real-rootedness
by scanning for complex zeros (min |B| off-axis / winding) over many (u1,u2,a,b)."""
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
def real_zeros(f,smax):
    zs=[]; s=mp.mpf('0.1'); prev=f(s); step=mp.mpf('0.15')
    while s<smax:
        s2=s+step; cur=f(s2)
        if prev*cur<0: zs.append(mp.findroot(f,(s+s2)/2))
        s,prev=s2,cur
    return zs
print("two-delta mixtures  B(s)=a·G((u1 s)²)+b·G((u2 s)²)  (a,b>0): complex zeros in strip?",flush=True)
print(" u1  u2   a   b | #real-zeros(s<15) | winding off-axis (0=LP, >real ⇒ COMPLEX ZEROS)",flush=True)
cases=[(1,2,1,1),(1,3,1,1),(1,5,1,1),(1,2,1,5),(1,2,5,1),(1,4,1,10),(1,1.5,1,1),(1,10,1,1),(1,3,3,1),(0.5,4,2,3)]
verdict_all_lp=True
for (u1,u2,a,b) in cases:
    u1=mp.mpf(u1);u2=mp.mpf(u2);a=mp.mpf(a);b=mp.mpf(b)
    f=lambda s: a*G((u1*s)**2)+b*G((u2*s)**2)
    rz=real_zeros(f,15)
    # winding over strip strictly above axis, up to Re=15, Im in [0.05,6]
    w=wind(f,mp.mpf('0.05'),15,mp.mpf('0.05'),6,600)
    lp = abs(w)<mp.mpf('0.5')
    verdict_all_lp=verdict_all_lp and lp
    print(f" {mp.nstr(u1,2)} {mp.nstr(u2,2)} {mp.nstr(a,2)} {mp.nstr(b,2)} | {len(rz):2d} real | wind={mp.nstr(w,3)} {'LP(no complex in strip)' if lp else '*** COMPLEX ZEROS — mixture-LP FALSE ***'}",flush=True)
print("",flush=True)
print(f"VERDICT: all two-delta mixtures LP in tested strip? {verdict_all_lp}",flush=True)
print("if TRUE across adversarial deltas ⇒ strong evidence mixture-LP holds ⇒ B∈LP UNCONDITIONAL (RH-free).",flush=True)
print("if FALSE ⇒ B's real-rootedness needs the modular Φ ⇒ the corner is NOT free, RH-adjacent.",flush=True)
print("DONE",flush=True)

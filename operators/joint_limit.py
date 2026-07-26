#!/usr/bin/env python3
"""§322: WHERE do zeta zeros enter the Jensen plane? Test the mollification-everywhere hypothesis:
bottom roots of J^{d,n} = zeros of B_n(s)=sum_k gamma_{n+k}(-s^2)^k/k! (shifted Borel transform, RH-free)
for ALL n — and zeta zeros appear at NO finite (n,d) window.
(1) B_0 zeros s_1..s_12; d=110 roots: sqrt(110)*t_i vs s_i vs gamma_i (any affine scaling) for i=1..12.
(2) n=6, n=16: B_n zeros vs sqrt(d)*t_i of J^{110,n}.
If Borel tracks everywhere and zeta nowhere: the mollification is uniform across the plane. Charisma: HANDLE."""
import mpmath as mp, pickle
mp.mp.dps=50
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam130.pkl","rb"))]
def Bn(n):
    def f(s):
        return sum(gam[n+k]*(-s*s)**k/mp.factorial(k) for k in range(min(100,130-n)))
    return f
def zeros_of(f,smax,count):
    zs=[]; s=mp.mpf('0.3'); prev=f(s); step=mp.mpf('0.2')
    while len(zs)<count and s<smax:
        s2=s+step; cur=f(s2)
        if prev*cur<0: zs.append(mp.findroot(f,(s+s2)/2))
        s,prev=s2,cur
    return zs
def jroots_t(d,n,take):
    co=[mp.binomial(d,k)*gam[n+k] for k in range(d,-1,-1)]
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=800,extraprec=400)
    return sorted(mp.sqrt(abs(mp.re(r))) for r in rts)[:take]
zeta=[mp.mpf(z) for z in ['14.134725','21.022040','25.010858','30.424876','32.935062','37.586178',
'40.918719','43.327073','48.005151','49.773832','52.970321','56.446248']]
d=110
for n in (0,6,16):
    print(f"--- n={n}, d={d} ---",flush=True)
    B=Bn(n)
    scount=12 if n==0 else 8
    szs=zeros_of(B,200,scount)
    ts=jroots_t(d,n,scount)
    print(" i : sqrt(d)*t_i    B_n zero s_i   ratio     gamma_i  (gamma ratios for contrast)",flush=True)
    for i in range(min(len(szs),len(ts))):
        lhs=mp.sqrt(d)*ts[i]
        print(f"{i+1:2d}: {mp.nstr(lhs,6):>10}  {mp.nstr(szs[i],6):>10}  {mp.nstr(lhs/szs[i],5):>7}   {mp.nstr(zeta[i],6):>8}  {mp.nstr(lhs/zeta[i],4)}",flush=True)
    # zeta-match test: best AFFINE map a*gamma_i+b fitted on i=1,2 -> check i>=3
    if len(ts)>=4:
        a=(mp.sqrt(d)*ts[1]-mp.sqrt(d)*ts[0])/(zeta[1]-zeta[0])
        b=mp.sqrt(d)*ts[0]-a*zeta[0]
        errs=[abs((a*zeta[i]+b)/(mp.sqrt(d)*ts[i])-1) for i in range(2,min(8,len(ts)))]
        print(f"  best-affine zeta fit residuals (i=3..): {[mp.nstr(e,3) for e in errs]}  (large => NOT zeta ladder)",flush=True)
    print("",flush=True)
print("DONE",flush=True)

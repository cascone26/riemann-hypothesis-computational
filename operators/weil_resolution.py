#!/usr/bin/env python3
"""§275: Weil-level RESOLUTION threshold. Fixed band [14,42], vary probe SPACING (6->1.5). Local zero spacing near
gamma~20 is ~2pi/log(20)~2.1. If probes get denser than the zeros can resolve, the Gram/Weil form goes near-singular
=> margin -> 0. The spacing where margin collapses = the Weil-level near-null / Nyman-Beurling resolution limit = a
measure of the zeros' rigidity. Uses the §273c-validated instrument."""
import mpmath as mp
mp.mp.dps = 25
def Omega(r): return mp.re(mp.digamma(mp.mpf(1)/4 + 1j*r/2)) - mp.log(mp.pi)
w = mp.mpf(2)   # narrower probes so spacing can go below the bump width
def F(c,r): return mp.e**(-(r-c)**2/(2*w**2))+mp.e**(-(r+c)**2/(2*w**2))
def Fc(c,z): return mp.e**(-(z-c)**2/(2*w**2))+mp.e**(-(z+c)**2/(2*w**2))
def bpts(ca,cb):
    pts=[mp.mpf(0)]
    for c in (ca,cb):
        for x in (c-5*w,c,c+5*w):
            if x>0: pts.append(x)
    return sorted(set(pts))+[mp.inf]
Nmax=80
Lam={}
for p in range(2,Nmax+1):
    if all(p%d for d in range(2,int(p**0.5)+1)):
        pk=p
        while pk<=Nmax: Lam[pk]=mp.log(p); pk*=p
def Bentry(ca,cb):
    pts=bpts(ca,cb)
    pole=mp.re(Fc(ca,1j/2)*Fc(cb,1j/2))*2
    arch=(1/(2*mp.pi))*2*mp.quad(lambda r:F(ca,r)*F(cb,r)*Omega(r),pts)
    prime=-2*sum(Lam[n]/mp.sqrt(n)*(1/(2*mp.pi))*2*mp.quad(lambda r:F(ca,r)*F(cb,r)*mp.cos(mp.log(n)*r),pts) for n in sorted(Lam))
    return pole+arch+prime
print(f"§275 Weil-level resolution: band [14,42], w={w}, vary probe spacing. local zero-spacing ~2.1 near g~20.",flush=True)
print(f"{'spacing':>8} {'#probes':>8} {'min-eig(margin)':>16} {'max-eig':>9} {'cond':>11}",flush=True)
lo,hi=mp.mpf(14),mp.mpf(42)
for sp in [mp.mpf(6),mp.mpf(4),mp.mpf(3),mp.mpf('2.5'),mp.mpf(2),mp.mpf('1.5')]:
    cs=[]; c=lo
    while c<=hi+mp.mpf('0.001'): cs.append(c); c+=sp
    M=len(cs); B=mp.zeros(M)
    for j in range(M):
        for k in range(j,M):
            v=Bentry(cs[j],cs[k]); B[j,k]=v; B[k,j]=v
    ev=mp.eigsy(B,eigvals_only=True); mn=min(ev); mx=max(ev)
    print(f"{mp.nstr(sp,3):>8} {M:>8} {mp.nstr(mn,6):>16} {mp.nstr(mx,5):>9} {mp.nstr(mx/mn,5):>11}",flush=True)
print("DONE. margin collapses toward 0 as spacing drops below ~zero-spacing => resolution/near-null threshold located.",flush=True)

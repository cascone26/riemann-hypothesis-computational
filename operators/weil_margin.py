#!/usr/bin/env python3
"""§274: Weil-positivity MARGIN study with the §273c-validated self-contained Weil instrument.
The full Weil form B (pole+arch-prime) = Sum_rho F_j(gamma)F_k(gamma) is PSD (RH). Its MIN eigenvalue = the positivity
margin. Sweep the number of zero-band centers (basis richness) and report min-eig: does the margin stay bounded, or
shrink toward 0 (a near-critical direction = Weil-level near-null / the zeros' resolution limit)?"""
import mpmath as mp
mp.mp.dps = 25
def Omega(r): return mp.re(mp.digamma(mp.mpf(1)/4 + 1j*r/2)) - mp.log(mp.pi)
w = mp.mpf(3)
allc = [mp.mpf(c) for c in [14,20,26,32,38,44,50,56,62,68]]
def F(c,r): return mp.e**(-(r-c)**2/(2*w**2))+mp.e**(-(r+c)**2/(2*w**2))
def Fc(c,z): return mp.e**(-(z-c)**2/(2*w**2))+mp.e**(-(z+c)**2/(2*w**2))
def bpts(ca,cb):
    pts=[mp.mpf(0)]
    for c in (ca,cb):
        for x in (c-4*w,c,c+4*w):
            if x>0: pts.append(x)
    return sorted(set(pts))+[mp.inf]
Nmax=60
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
print(f"§274 Weil-positivity margin vs basis richness (validated instrument, w={w}, primes<={Nmax}, dps={mp.mp.dps})",flush=True)
print(f"{'#centers':>9} {'band':>10} {'min-eig(margin)':>16} {'max-eig':>10} {'cond':>10}",flush=True)
for M in range(2, len(allc)+1):
    cs=allc[:M]
    B=mp.zeros(M)
    for j in range(M):
        for k in range(j,M):
            v=Bentry(cs[j],cs[k]); B[j,k]=v; B[k,j]=v
    ev=mp.eigsy(B,eigvals_only=True)
    mn=min(ev); mx=max(ev)
    print(f"{M:>9} {'['+str(int(cs[0]))+','+str(int(cs[-1]))+']':>10} {mp.nstr(mn,6):>16} {mp.nstr(mx,5):>10} {mp.nstr(mx/mn,4):>10}",flush=True)
print("DONE. margin bounded away from 0 => robust Weil positivity; margin -> 0 => near-critical direction (Weil near-null).",flush=True)

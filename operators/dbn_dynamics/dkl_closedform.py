#!/usr/bin/env python3
"""§199.3b step1: validate the hand-derived closed form of D_kl(t)=C_kl+C_lk (cos-overlap, t∈[0,2a]) vs numeric.
β=π/2a. k,l odd. Derived (product-to-sum, using 2(k±l)βa=(k±l)π and k±l even):
  k≠l: D_kl = klβ{ [sin(lβt)-sin(kβt)]/(k-l) - [sin(lβt)+sin(kβt)]/(k+l) }
  k=l: D_kk = k²β²(2a-t)cos(kβt) - kβ sin(kβt)
Validate vs numeric C_kl+C_lk (nested quad)."""
import mpmath as mp
mp.mp.dps=30
a=mp.mpf('0.549'); b=mp.pi/(2*a)
def dprime(k,x): return (k*mp.pi/(2*a))*mp.cos(k*mp.pi*(x+a)/(2*a))
def Cnum(k,l,t): return mp.quad(lambda x: dprime(k,x)*dprime(l,x-t),[t-a,a])
def Dnum(k,l,t): return Cnum(k,l,t)+Cnum(l,k,t)
def Dcf(k,l,t):
    if k==l: return k*k*b*b*(2*a-t)*mp.cos(k*b*t) - k*b*mp.sin(k*b*t)
    return k*l*b*( (mp.sin(l*b*t)-mp.sin(k*b*t))/(k-l) - (mp.sin(l*b*t)+mp.sin(k*b*t))/(k+l) )
print("§199.3b.1 D_kl closed-form vs numeric  @ a=0.549")
maxe=mp.mpf(0)
for (k,l) in [(1,1),(1,3),(3,3),(1,5),(3,7),(5,7),(7,7)]:
    for t in [mp.mpf('0.2'),mp.mpf('0.5'),mp.mpf('0.9')]:
        dn=Dnum(k,l,t); dc=Dcf(k,l,t); e=abs(dn-dc); maxe=max(maxe,e)
        print(f"  D[{k},{l}]({mp.nstr(t,2)}): num={mp.nstr(dn,10)} cf={mp.nstr(dc,10)} err={mp.nstr(e,3)}")
print(f"\nmax err = {mp.nstr(maxe,3)}  (≈0 ⟹ D_kl closed form CORRECT — safe to build entry integrals on it)")

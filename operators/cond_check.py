import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
import mpmath as mp
mp.mp.dps=18
from flint import ctx
import dh_certified_tau as DH
def psiR(x,L,nu,cond):
    x=mp.mpf(x);nu=mp.mpf(nu);L=mp.mpf(L);logc=mp.log(mp.pi/cond)
    hp=lambda r: mp.re(mp.digamma(nu+1j*r/2))-logc
    def Sr(r):
        B=2*mp.pi*x/L;den=B*B-r*r
        if abs(den)<mp.mpf(10)**-18: r=r+mp.mpf(10)**-12;den=B*B-r*r
        return (mp.cos(r*L)-mp.cos(2*mp.pi*x))*B/den
    f=lambda r: hp(r)*Sr(r);B=abs(2*mp.pi*x/L);per=mp.pi/L;R=mp.mpf(50);pts=[mp.mpf(0)];k=1
    while k*per<R: pts.append(k*per);k+=1
    pts.append(R);allp=sorted(set([float(p) for p in pts]+[float(B)]))
    val=2*mp.quad(f,allp)+2*B*mp.cos(2*mp.pi*x)*mp.quad(lambda r: hp(r)/r**2,[R,4*R,mp.inf])
    return -val/(2*mp.pi**2)
c=14;prec=400;ctx.prec=prec
S,CC,XC,L=DH.closed_forms_nu(6,c,prec,0.25); Lf=float(L.mid());Lm=float(mp.log(c));h=mp.mpf(10)**-5
def kj(nu,cond):
    n=1; d=float((psiR(n+h,Lm,nu,cond)-psiR(n-h,Lm,nu,cond))/(2*h))
    return d-2*float(CC[n].mid())+(2/Lf)*float(XC[n].mid())
import math
kj_c1=kj(mp.mpf(1)/4,1); kj_c5=kj(mp.mpf(1)/4,5)
print(f"(kappa+J)[nu=1/4,cond=1] numeric = {kj_c1:.5f}  (Groskin exact 4.30204)")
print(f"(kappa+J)[nu=1/4,cond=5] numeric = {kj_c5:.5f}")
print(f"difference = {kj_c5-kj_c1:.5f}   log(5) = {math.log(5):.5f}   => conductor term -log(pi/cond) {'CONFIRMED' if abs((kj_c5-kj_c1)-math.log(5))<0.01 else 'NO'}")

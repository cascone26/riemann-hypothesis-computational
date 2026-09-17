import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
import mpmath as mp
mp.mp.dps=20
from flint import ctx
import dh_certified_tau as DH
from dh_certified_tau import build_tau_general, build_arb_tau, margin, dh_prime_side

def psiR(x,L,nu,cond):
    x=mp.mpf(x);nu=mp.mpf(nu);L=mp.mpf(L);logc=mp.log(mp.pi/cond)
    hp=lambda r: mp.re(mp.digamma(nu+1j*r/2))-logc
    def Sr(r):
        B=2*mp.pi*x/L;den=B*B-r*r
        if abs(den)<mp.mpf(10)**-20: r=r+mp.mpf(10)**-14;den=B*B-r*r
        return (mp.cos(r*L)-mp.cos(2*mp.pi*x))*B/den
    f=lambda r: hp(r)*Sr(r);B=abs(2*mp.pi*x/L)
    per=mp.pi/L;R=mp.mpf(60);pts=[mp.mpf(0)];k=1
    while k*per<R: pts.append(k*per);k+=1
    pts.append(R);allp=sorted(set([float(p) for p in pts]+[float(B)]))
    val=2*mp.quad(f,allp)+2*B*mp.cos(2*mp.pi*x)*mp.quad(lambda r: hp(r)/r**2,[R,4*R,mp.inf])
    return -val/(2*mp.pi**2)

def kj_dh(c,prec):
    S,CC,XC,L=DH.closed_forms_nu(8,c,prec,0.75); Lf=float(L.mid());Lm=float(mp.log(c));h=mp.mpf(10)**-5
    vals=[]
    for n in [1,2,3]:
        d=float((psiR(n+h,Lm,mp.mpf(3)/4,5)-psiR(n-h,Lm,mp.mpf(3)/4,5))/(2*h))
        vals.append(d-2*float(CC[n].mid())+(2/Lf)*float(XC[n].mid()))
    return sum(vals)/3, min(vals), max(vals)

prec,N=800,64
print(f"MONOTONE TABLE: certified D-H n_neg vs zeta n_neg, growing prime cutoff c (N={N})",flush=True)
print(f"{'c':>4} | {'(k+J)_DH':>18} | {'zeta n_neg':>10} | {'D-H n_neg [at kj, kj-.005, kj+.005]':>36}",flush=True)
for c in [25,75,100]:
    ctx.prec=prec
    kjc,klo,khi=kj_dh(c,prec)
    dw,dpos=dh_prime_side(c,prec)
    Az,Dz=build_arb_tau(c,N,prec); nz,_,_=margin(Az,Dz)
    ns=[]
    for kj in [kjc,kjc-0.005,kjc+0.005]:
        A,DIM=build_tau_general(c,N,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=dw,prime_positions=dpos)
        nn,_,_=margin(A,DIM); ns.append(nn)
    print(f"{c:>4} | {kjc:.5f} [{klo:.4f},{khi:.4f}] | {nz:>10} | {ns}",flush=True)

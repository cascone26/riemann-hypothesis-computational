import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
import mpmath as mp
mp.mp.dps=20
from flint import arb, ctx
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

c=50;prec=800;ctx.prec=prec
S,CC,XC,L=DH.closed_forms_nu(8,c,prec,0.75); Lf=float(L.mid());Lm=float(mp.log(c));h=mp.mpf(10)**-5
KJs=[]
for n in [1,2]:
    d=float((psiR(n+h,Lm,mp.mpf(3)/4,5)-psiR(n-h,Lm,mp.mpf(3)/4,5))/(2*h))
    kj=d-2*float(CC[n].mid())+(2/Lf)*float(XC[n].mid()); KJs.append(kj)
    print(f"  c=50 n={n}: (kappa+J)_DH={kj:.6f}",flush=True)
kjc=sum(KJs)/len(KJs)
print(f"  => (kappa+J)_DH(c=50) ~ {kjc:.5f}",flush=True)
dw,dpos=dh_prime_side(c,prec)
print(f"D-H vs zeta at c=50, N=64 (sweep constant +/-0.01 for robustness):",flush=True)
Az,Dz=build_arb_tau(c,64,prec); nz,_,_=margin(Az,Dz)
print(f"  zeta(exact): n_neg={nz}",flush=True)
for dk in [-0.01,-0.005,0.0,0.005,0.01]:
    kj=kjc+dk
    A,DIM=build_tau_general(c,64,prec,nu=0.75,include_pole=False,kj_const=kj,prime_weights=dw,prime_positions=dpos)
    nn,pd,mg=margin(A,DIM)
    print(f"  D-H (k+J={kj:.4f}): n_neg={nn}",flush=True)

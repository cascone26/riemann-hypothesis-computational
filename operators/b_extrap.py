"""Pin the constant b in lambda_a = c log(1/a) + b + O(a) via fast (numpy) archimedean lambda_a
at many small a + Richardson-style extrapolation (a->0 kills O(a)). Then closed-form check."""
import numpy as np, os
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,"suzuki_screw_lambda_a.py")).read()
cut=src.find("if __name__"); exec(src[:cut] if cut>0 else src, globals())
# archimedean-only lambda_a for small a (primes absent for a<0.347); use screw_g_no_primes if present
def lam_arch(a, M=6000, K=50):
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    try:
        gv=np.array([screw_g_no_primes(d, None, 0) for d in diffs])
    except Exception:
        gv=np.array([screw_g(d, vonmangoldt_upto(2), 2) for d in diffs])
    idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Q=Wp.T@G@Wp; return (np.linalg.eigvalsh((Q+Q.T)/2)/a)[0]
avs=np.array([0.005,0.008,0.012,0.018,0.026,0.038,0.055])
las=np.array([lam_arch(a) for a in avs])
for a,l in zip(avs,las): print(f"a={a:.3f}: lam={l:.6f}  lam-log(1/a)={l-np.log(1/a):.6f}")
# fit c log(1/a)+b+e*a using smallest a's
X=np.column_stack([np.log(1/avs),np.ones_like(avs),avs])
c,b,e=np.linalg.lstsq(X,las,rcond=None)[0]
print(f"\nfit: c={c:.5f} (Suzuki=1)  b={b:.6f}  e={e:.4f}")
# also 2-param on 3 smallest (less O(a) leverage)
X2=np.column_stack([np.log(1/avs[:3]),np.ones(3)])
c2,b2=np.linalg.lstsq(X2,las[:3],rcond=None)[0]
print(f"2-param (3 smallest): c={c2:.4f} b={b2:.5f}")
const=0.5*(np.log(2*np.pi)-(-np.euler_gamma-2*np.log(2)+2))  # psi(2)=1-gamma... = -gamma+1
# psi(2) = 1 - gamma
psi2=1-np.euler_gamma
const=0.5*(np.log(2*np.pi)-psi2)
print(f"\nb={b:.5f}; (log2pi-psi2)/2={const:.5f}; mu1=b-const={b-const:.5f}; b+const={b+const:.5f}")
print("closed-form checks for b:")
import math
for nm,v in [("-2",-2),("-(1+gamma)",-(1+np.euler_gamma)),("-log(2pi)",-np.log(2*np.pi)),
             ("-(gamma+log2)",-(np.euler_gamma+np.log(2))),("psi2-2",psi2-2),
             ("-2+gamma/2",-2+np.euler_gamma/2),("-(2-1/4)",-(2-0.25))]:
    print(f"  b vs {nm}={v:.5f}: diff={abs(b-v):.5f}")
print("DONE")

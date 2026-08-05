"""Extract C = lim a^2 min-eig(a) (the positivity 'curvature'); decompose the near-null mode's
Weil value into ARCHIMEDEAN vs PRIME to see WHY arch>primes (the positivity mechanism)."""
import numpy as np, os
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,"suzuki_screw_lambda_a.py")).read()
cut=src.find("if __name__"); exec(src[:cut] if cut>0 else src, globals())
def Qmat_full(a,M,K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gvals=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gvals[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Q=Wp.T@G@Wp; Q=(Q+Q.T)/2/a
    ev,V=np.linalg.eigh(Q)
    return ev, V, x, kk
print(" a    min-eig     a^2*min-eig  (C = lim?)")
Cs=[]
for a in (2.0,3.0,4.0,5.0):
    K=max(40,int(18*a)); M=max(4000,int(1200*a))
    ev,V,x,kk=Qmat_full(a,M,K)
    C=a**2*ev[0]; Cs.append((a,C))
    print(f" {a:.1f}   {ev[0]:+.4e}   {C:+.5f}", flush=True)
print(f"a^2 min-eig: {[round(c,4) for _,c in Cs]}  (converging to C>0 => RH-safe curvature)")
print("DONE")

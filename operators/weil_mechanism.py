"""Positivity MECHANISM: decompose the near-null (smoothest) mode's Weil value into
ARCHIMEDEAN vs PRIME. W(phi0)=arch-prime>0 is why RH holds; track how arch>prime as a grows."""
import numpy as np, os
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,"suzuki_screw_lambda_a.py")).read()
cut=src.find("if __name__"); exec(src[:cut] if cut>0 else src, globals())
# screw_g(d) = sum_n Lambda(n)-kernel (the PRIME part of the kernel). The Weil form Q built from it
# is the FULL form? Check: is screw_g the prime kernel or the full (arch-prime)? Inspect its zero-arg.
Lam=vonmangoldt_upto(1000)
print("screw_g(0) =", screw_g(0.0, Lam, 1000), "(kernel diagonal)")
print("screw_g(log2) =", screw_g(np.log(2), Lam, 1000))
# The near-null mode at a=4: get it, and compute its Rayleigh quotient decomposition
def Qmat_V(a,M,K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    L=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,L,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Q=Wp.T@G@Wp; Q=(Q+Q.T)/2/a
    ev,V=np.linalg.eigh(Q)
    return ev,V,x,kk,Wp,G
for a in (3.0,4.0,5.0):
    ev,V,x,kk,Wp,G=Qmat_V(a,int(1200*a),int(18*a))
    v0=V[:,0]
    # near-null mode value and its localization (is it smoothest = dominated by lowest kk?)
    dom=np.argmax(np.abs(v0))+1
    print(f"a={a}: min-eig={ev[0]:+.4e}, near-null dominant mode k={dom}/{len(kk)}, "
          f"low-freq weight |v0[0:3]|={np.linalg.norm(v0[:3])/np.linalg.norm(v0):.3f}", flush=True)
print("DONE")

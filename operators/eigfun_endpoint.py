import numpy as np
from suzuki_screw_lambda_a import screw_g_no_primes
import suzuki_screw_lambda_a as S
S.screw_g=screw_g_no_primes
# reconstruct arch eigenfunction v0 at a, on a grid, from sine-basis coefficients
def eigfun(a,M,K):
    from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in diffs]); idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg); w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]; Q=Wp.T@G@Wp; Q=(Q+Q.T)/2
    val,vec=np.linalg.eigh(Q); c=vec[:,0]
    V=np.sin(arg)@c   # v0(x) = sum c_k sin(kπ(x+a)/2a)
    return x,V
a=0.36; x,V=eigfun(a,7000,120)
V=V/np.max(np.abs(V))
# test v0 ~ (a^2-x^2)^{1/2} near endpoint: fit log|V| vs log(a-x) for x near a
mask=(x>0.85*a)&(x<0.999*a)
p=np.polyfit(np.log(a-x[mask]),np.log(np.abs(V[mask])),1)
print(f"eigenfunction endpoint fit (a={a}): |v0(x)| ~ (a-x)^{p[0]:.3f} near x=a",flush=True)
print(f"  ⟹ v0 ~ (a-x)^{p[0]:.3f}, so v0' ~ (a-x)^{p[0]-1:.3f}  [predict exponent 1/2 ⟹ v0'~(a-x)^-1/2, beta=1/2]",flush=True)
# also compare to sqrt directly
pred=np.sqrt(a**2-x[mask]**2); pred=pred/np.max(np.abs(pred))*np.max(np.abs(V[mask]))
print(f"  correlation of v0 with (a^2-x^2)^1/2 near endpoint: {np.corrcoef(np.abs(V[mask]),np.abs(pred))[0,1]:.5f}",flush=True)

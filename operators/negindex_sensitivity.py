import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
def boundary_kernel(t): return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)
def build_Q(a,M,K,kernel):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([kernel(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    G=gf[shift]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]; Q=Wp.T@G@Wp; return (Q+Q.T)/2
a,M,K=1.0,3000,80
# 1) real screw form
Qf=build_Q(a,M,K,screw_g); nf=int((np.linalg.eigvalsh(Qf)<-1e-9).sum())
# 2) inject a synthetic OFF-LINE zero: add its explicit-formula contribution to the kernel.
#    An off-line zero rho=1/2+delta+i*gamma contributes -2 cos(gamma t) e^{-delta|t|}... to the density;
#    as a screw/Weil kernel perturbation on g we add a term that BREAKS positivity if delta!=0.
def screw_offline(t,Lam,Nmax,delta=0.3,gamma=14.13):
    base=screw_g(t,Lam,Nmax)
    # a pair {rho,1-bar rho} off-line adds a real oscillatory growing bump ~ +/- e^{delta|t|} cos(gamma t)
    return base + 4.0*np.cosh(delta*abs(t))*np.cos(gamma*t) - 4.0*np.cos(gamma*t)
Qo=build_Q(a,M,K,lambda t,L,N: screw_offline(t,L,N)); no=int((np.linalg.eigvalsh(Qo)<-1e-9).sum())
print(f"a={a}: #neg(real screw Q_full) = {nf}   (arsenal: 0)")
print(f"a={a}: #neg(off-line-injected Q)= {no}   (if >0, our #neg IS RH-sensitive, NOT AF-vacuous)")
ev=np.linalg.eigvalsh(Qo); print(f"   off-line lowest 4 eigs: {np.round(ev[:4],4)}")

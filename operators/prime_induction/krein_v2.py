# Krein extension, vectorized: objective = min S over xi>=2 (0-well handled separately).
import numpy as np, sys, os
from scipy.optimize import minimize
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
from scipy.special import digamma
def run(a, T=3.0, ncp=5):
    xis=np.linspace(2.0,150,3500)
    Om=np.real(digamma(0.25+0.5j*xis))-np.log(np.pi)
    ts1=np.linspace(0,2*a,2500); dt1=ts1[1]-ts1[0]; w1=np.full_like(ts1,dt1); w1[0]=w1[-1]=dt1/2
    ts2=np.linspace(2*a,2*a+T,2500); dt2=ts2[1]-ts2[0]; w2=np.full_like(ts2,dt2); w2[0]=w2[-1]=dt2/2
    k1=-4*(np.exp(ts1/2)+np.exp(-ts1/2)-2); k2a=-4*(np.exp(a)+np.exp(-a)-2)
    C1=np.cos(np.outer(xis,ts1)); C2=np.cos(np.outer(xis,ts2))
    base=Om+2*(xis**2)*(C1@(w1*k1))
    Nmax=int(np.floor(np.exp(2*a)))
    if Nmax>=2:
        Lam=S.vonmangoldt_upto(max(Nmax,3))
        for n in range(2,Nmax+1):
            if Lam[n]>0: base += -2*(Lam[n]/np.sqrt(n))*np.cos(xis*np.log(n))
    uu2=(ts2-2*a)/T; ucp=np.linspace(0,1,ncp+2)
    def Sof(p):
        cp=np.concatenate([[1.0],p,[0.0]]); f=np.interp(uu2,ucp,cp)
        return base+2*(xis**2)*(C2@(w2*(k2a*f)))
    def neg(p): return -np.min(Sof(p))
    best=None
    for seed in range(6):
        rng=np.random.RandomState(seed)
        p0=np.clip(np.linspace(0.85,0.1,ncp)+0.15*rng.randn(ncp),-0.2,1.4)
        r=minimize(neg,p0,method='Nelder-Mead',options={'maxiter':1500,'xatol':1e-4,'fatol':1e-6})
        if best is None or r.fun<best.fun: best=r
    Sv=Sof(best.x); i=np.argmin(Sv)
    print(f"a={a:.2f}: min_(xi>=2) S = {-best.fun:+.4f} at xi={xis[i]:.1f}  cp={np.round(best.x,2)}", flush=True)
for a in [0.30,0.37,0.42,0.48]: run(a)
print("\nmin>=0 => extension certificate exists off the 0-well at that a (=> proof shape: certificate + explicit 0-well estimate).")

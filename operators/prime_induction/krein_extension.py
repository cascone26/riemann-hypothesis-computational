# KREIN EXTENSION route: find an extension of the pole kernel beyond |t|=2a making the TOTAL multiplier
# S_a(xi) >= 0 everywhere. If found => lambda_a >= 0 PROVEN directly (Q=(1/2pi)int S|vhat|^2).
# Parametrize the pole-kernel tail on [2a, 2a+T] as k(2a)*sum spline/decay basis; optimize to maximize min_xi S.
import numpy as np, sys, os
from scipy.optimize import minimize
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
from scipy.special import digamma
def Omega(xi): return np.real(digamma(0.25+0.5j*xi))-np.log(np.pi)
def total_S(a, params, xis, T=3.0, nb=6):
    # pole kernel exact on [0,2a]; tail on [2a,2a+T] = k(2a) * f(u), f(0)=1, f(T)=0, f = interp of params (>=0 free)
    ts1=np.linspace(0,2*a,3000); ts2=np.linspace(2*a,2*a+T,3000)
    k1=-4*(np.exp(ts1/2)+np.exp(-ts1/2)-2)
    k2a=-4*(np.exp(a)+np.exp(-a)-2)
    # tail shape: cubic-interpolated control points from 1 -> 0
    cp=np.concatenate([[1.0],params,[0.0]])
    uu=np.linspace(0,1,len(cp))
    f=np.interp((ts2-2*a)/T,uu,cp)
    k2=k2a*f
    def cosFT(ts,k,xi):
        dt=ts[1]-ts[0]; w=np.full_like(ts,dt); w[0]=w[-1]=dt/2
        return np.sum(w*k*np.cos(xi*ts))
    out=np.empty_like(xis)
    for i,xi in enumerate(xis):
        out[i]=Omega(xi)+2*(xi**2)*(cosFT(ts1,k1,xi)+cosFT(ts2,k2,xi))
    # primes
    Nmax=int(np.floor(np.exp(2*a)))
    if Nmax>=2:
        Lam=S.vonmangoldt_upto(max(Nmax,3))
        for n in range(2,Nmax+1):
            if Lam[n]>0: out += -2*(Lam[n]/np.sqrt(n))*np.cos(xis*np.log(n))
    return out
for a in [0.30,0.37,0.42]:
    xis=np.linspace(0,120,3000)
    def neg_minS(p): return -np.min(total_S(a,p,xis))
    best=None
    for seed in range(3):
        rng=np.random.RandomState(seed)
        p0=np.clip(np.linspace(0.8,0.1,4)+0.1*rng.randn(4),0,1.2)
        r=minimize(neg_minS,p0,method='Nelder-Mead',options={'maxiter':400,'xatol':1e-3,'fatol':1e-5})
        if best is None or r.fun<best.fun: best=r
    minS=-best.fun
    print(f"a={a:.2f}: optimized extension -> min_xi S_a = {minS:+.4f}   {'>=0 => LAMBDA_a>=0 PROVEN (numerically)' if minS>=0 else '(still dips)'}", flush=True)
print("\nS>=0 achieved => Krein-extension certificate exists at that a => first-step positivity via explicit certificate.")

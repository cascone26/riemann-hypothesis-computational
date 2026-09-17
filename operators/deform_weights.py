import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from dh_certified_tau import build_tau_general, margin
from _dh_negative_control import lambda_D

def vonmangoldt(N):
    lam=[0.0]*(N+1); comp=[False]*(N+1)
    for p in range(2,N+1):
        if not comp[p]:
            for m in range(2*p,N+1,p): comp[m]=True
            lp=math.log(p);pk=p
            while pk<=N: lam[pk]=lp;pk*=p
    return lam

c,N,prec=50,64,800
ctx.prec=prec
VM=vonmangoldt(c); LD=lambda_D(c)
pos=[arb(n).log() for n in range(2,c+1)]
print(f"Fix zeta arch (nu=1/4,pole); deform prime weights w(t)=(1-t)Lambda + t Lambda_D, all n<=c={c}, N={N}")
print(f"{'t':>6} | {'n_neg':>6}")
for t in [0.0,0.02,0.05,0.1,0.2,0.35,0.5,0.75,1.0]:
    w=[arb((1-t)*VM[n] + t*float(LD[n]))*(arb(n)**arb('-0.5')) for n in range(2,c+1)]
    A,D=build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=w,prime_positions=pos)
    nn,_,_=margin(A,D)
    print(f"{t:>6.2f} | {nn:>6}",flush=True)
print("t=0 must be n_neg=0 (zeta control). First t with n_neg>0 = positivity boundary in the non-mult direction.")

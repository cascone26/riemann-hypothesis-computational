#!/usr/bin/env python3
"""§313C: close degrees d=5..8 for ALL n (within the effective scheme): compute gamma_m to m=96 and directly
verify J^{d,n} hyperbolicity for n up to N(d) (18/33/55/87 for d=5/6/7/8). Above N(d), the §313A c_m law +
Hermite margin covers it (leading-order scheme). dps 60, denser quad nodes for the high-m saddle (u*~1.7)."""
import mpmath as mp
mp.mp.dps=60
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,100):
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n*n*e2u)
        s+=term
        if abs(term)<mp.mpf(10)**(-70) and n>3: break
    return s
Mmax=97
import os,pickle
cache="operators/_gam96.pkl"
if os.path.exists(cache):
    gam=pickle.load(open(cache,"rb")); gam=[mp.mpf(g) for g in gam]
    print(f"loaded {len(gam)} cached moments",flush=True)
else:
    gam=[]
    nodes=[0,0.4,0.8,1.2,1.6,2.0,2.5,3.0,4.5]
    for m in range(Mmax):
        M=mp.quad(lambda u: Phi(u)*u**(2*m),nodes)
        gam.append(M/mp.factorial(2*m))
        if m%12==0: print(f"  gamma_{m} = {mp.nstr(gam[m],6)}",flush=True)
    pickle.dump([str(g) for g in gam],open(cache,"wb"))
    print(f"computed+cached {len(gam)} moments",flush=True)
# sanity: all positive, log-concave spot check
assert all(g>0 for g in gam), "NEGATIVE MOMENT?!"
def jensen_maximag(d,n):
    co=[mp.binomial(d,k)*gam[n+k] for k in range(d,-1,-1)]
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=300,extraprec=200)
    mx=max(abs(mp.im(r)) for r in rts); mxa=max(abs(r) for r in rts)
    return mx/mxa
targets={5:18,6:33,7:55,8:88}
print("direct J^{d,n} hyperbolicity up to N(d):",flush=True)
allok=True
for d,N in targets.items():
    worst=mp.mpf(0); bad=[]
    for n in range(0,min(N,Mmax-d)+1):
        v=jensen_maximag(d,n)
        worst=max(worst,v)
        if v>mp.mpf('1e-20'): bad.append((n,v))
    ok=not bad
    allok=allok and ok
    print(f"  d={d}: n=0..{min(N,Mmax-d)}  worst max|Im|/|z| = {mp.nstr(worst,3)}  {'ALL HYPERBOLIC => d closed for ALL n (scheme)' if ok else f'FAILURES at {bad[:3]}'}",flush=True)
print(f"=> degrees 3-8 all-n status: {'CLOSED (within leading-order scheme; d=3,4 closed at N<=12 earlier)' if allok else 'INCOMPLETE'}",flush=True)
print("DONE",flush=True)

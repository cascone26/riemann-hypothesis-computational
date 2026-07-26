#!/usr/bin/env python3
"""§314b: finish d=8 under the proved margin — extend moments to m=106, check J^{8,n} n=89..97."""
import mpmath as mp, pickle, os
mp.mp.dps=60
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,100):
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n*n*e2u)
        s+=term
        if abs(term)<mp.mpf(10)**(-70) and n>3: break
    return s
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam96.pkl","rb"))]
nodes=[0,0.4,0.8,1.2,1.6,2.0,2.5,3.0,4.5]
for m in range(len(gam),107):
    M=mp.quad(lambda u: Phi(u)*u**(2*m),nodes)
    gam.append(M/mp.factorial(2*m))
    assert gam[-1]>0
print(f"moments extended to m={len(gam)-1} (all >0)",flush=True)
pickle.dump([str(g) for g in gam],open("operators/_gam106.pkl","wb"))
def jensen_maximag(d,n):
    co=[mp.binomial(d,k)*gam[n+k] for k in range(d,-1,-1)]
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=300,extraprec=200)
    return max(abs(mp.im(r)) for r in rts)/max(abs(r) for r in rts)
worst=mp.mpf(0)
for n in range(89,98):
    worst=max(worst,jensen_maximag(8,n))
print(f"d=8, n=89..97: worst max|Im|/|z| = {mp.nstr(worst,3)}  => {'d=8 CLOSED for ALL n (proved-margin scheme)' if worst<mp.mpf('1e-20') else 'FAILURE'}",flush=True)
print("DONE",flush=True)

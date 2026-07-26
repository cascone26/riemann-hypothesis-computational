#!/usr/bin/env python3
"""§316: directional ladder through d=25. For d=9..25: coarse tau_d^+ (dps35), N_dir(d)=1/sqrt(3 tau_d),
extend moment cache to m=130, verify J^{d,n} hyperbolic for n=0..N_dir(d). Closes each d for ALL n under the
directional scheme (leading-order; rigor tier as §315). Charisma: computations HANDLE, scheme REPORT."""
import mpmath as mp, pickle, os
mp.mp.dps=35
def hermite_coeffs(d):
    P=[mp.mpf(1)]; Q=[mp.mpf(0),mp.mpf(1)]
    for k in range(1,d):
        R=[mp.mpf(0)]+Q
        for i in range(len(P)): R[i]-=k*P[i]
        P,Q=Q,R
    return Q
def is_hyp(co,tol='1e-12'):
    try: rts=mp.polyroots(list(reversed(co)),maxsteps=300,extraprec=160)
    except Exception: return False
    return max(abs(mp.im(r)) for r in rts)/max(abs(r) for r in rts) < mp.mpf(tol)
def tau_plus(d):
    co=hermite_coeffs(d); h=mp.mpf(d)/2
    g=lambda k:(k-h)**3-(3*h*h/4)*(k-h)
    hyp=lambda t: is_hyp([c*mp.e**(t*g(k)) for k,c in enumerate(co)])
    lo,hi=mp.mpf(0),mp.mpf('0.05')
    while hyp(hi) and hi<64: hi*=2
    if hi>=64: return mp.mpf(64)
    for _ in range(26):
        mid=(lo+hi)/2
        if hyp(mid): lo=mid
        else: hi=mid
    return lo
taus={}; Nd={}
print("d : tau_d^+       N_dir(d)",flush=True)
for d in range(9,26):
    taus[d]=tau_plus(d)
    Nd[d]=int(mp.ceil(1/mp.sqrt(3*taus[d])))
    print(f"{d:2d}: {mp.nstr(taus[d],5):>10}  {Nd[d]}",flush=True)
maxm=max(Nd[d]+d for d in Nd)
print(f"max moment needed: {maxm}",flush=True)
# extend moments
mp.mp.dps=60
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,100):
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n*n*e2u)
        s+=term
        if abs(term)<mp.mpf(10)**(-70) and n>3: break
    return s
src="operators/_gam130.pkl" if os.path.exists("operators/_gam130.pkl") else "operators/_gam106.pkl"
gam=[mp.mpf(g) for g in pickle.load(open(src,"rb"))]
nodes=[0,0.4,0.8,1.2,1.6,2.0,2.5,3.0,4.5]
target=max(maxm,130)
for m in range(len(gam),target+1):
    M=mp.quad(lambda u: Phi(u)*u**(2*m),nodes)
    gam.append(M/mp.factorial(2*m))
    assert gam[-1]>0
pickle.dump([str(g) for g in gam],open("operators/_gam130.pkl","wb"))
print(f"moments to m={len(gam)-1} (all >0, cached)",flush=True)
def jensen_maximag(d,n):
    co=[mp.binomial(d,k)*gam[n+k] for k in range(d,-1,-1)]
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=400,extraprec=250)
    return max(abs(mp.im(r)) for r in rts)/max(abs(r) for r in rts)
print("direct checks below N_dir(d):",flush=True)
allok=True
for d in range(9,26):
    worst=mp.mpf(0)
    for n in range(0,Nd[d]+1):
        worst=max(worst,jensen_maximag(d,n))
    ok=worst<mp.mpf('1e-18')
    allok=allok and ok
    print(f"  d={d}: n=0..{Nd[d]}  worst|Im|/|z|={mp.nstr(worst,3)}  {'CLOSED all n (directional scheme)' if ok else 'FAIL?!'}",flush=True)
print(f"=> degrees 9-25: {'ALL CLOSED under directional ladder' if allok else 'INCOMPLETE'}",flush=True)
print("DONE",flush=True)

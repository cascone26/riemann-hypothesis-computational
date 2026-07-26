#!/usr/bin/env python3
"""§326: the d*(T0) EXCHANGE-RATE EXPERIMENT (Roundtable-converged, sealed protocol).
PRE-REGISTERED (fixed before any fit):
 FAMILY: Xi_T0(z) = Xi_xi(z) * (1 - z^2/a^2)(1 - z^2/abar^2), a = T0 + i*DELTA, ONE off-line pair
   planted at height T0, fixed depth DELTA. Xi_xi = the true zeta xi (cached moments _gam130.pkl, hyperbolic
   ladder verified to d>=25 + published to ~10^17). KNOB: T0. FIXED: DELTA=1.0.
 MOMENTS of the target (derivation in header): with beta=-2Re(1/a^2), kappa=|a|^-4,
   a_m = gamma_m - beta*(2m)(2m-1)*gamma_{m-1} + kappa*(2m)(2m-1)(2m-2)(2m-3)*gamma_{m-2}
   (beta<0 for T0>DELTA => all-positive sequence). a_m are the Jensen sequence (matches lp_ladder25 convention).
 MEASUREMENT: d*(T0) = smallest degree d s.t. some Jensen poly J^{d,n}(X)=sum_k C(d,k) a_{n+k} X^k
   (n=0..Mmax-d) is NON-hyperbolic (has a complex root). TWO-SIDED cert (Euler): report the failing (d*,n*)
   with its max|Im|/|z| ABOVE threshold AND confirm d*-1 fully hyperbolic (all n) BELOW threshold.
 NULL (pre-registered): d*(T0) = c * T0^2  (GORTTW Thm 1.2 resolution rate). Fit log d* vs log T0; report
   exponent WITH residuals across the T0 span (Gauss: 2 points is numerology).
 INVARIANCE STAMP: Jensen margins invariant under z->lam z & Xi->c Xi; base scale fixed by the real xi, so T0
   is in xi's own height units (dimensionless against the zeta zero scale gamma_1~14.13). DELTA=1 is a fixed
   fraction of that scale. (No hidden gauge freedom - the base is a fixed real function, not a rescalable Gaussian.)
 MEMBERSHIP verdict (logged separately, NOT a gate): Xi_T0 is off-line by construction (that's the point);
   whether it comes from a positive kernel is irrelevant to the Jensen-ladder measurement.
"""
import mpmath as mp, pickle, math
mp.mp.dps=40
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam130.pkl","rb"))]
Mmax=len(gam)  # 131
DELTA=mp.mpf(1)
THR=mp.mpf('1e-11')   # hyperbolicity threshold on max|Im z|/max|z|
def target_moments(T0):
    T0=mp.mpf(T0)
    a2r=T0*T0-DELTA*DELTA; a2i=2*T0*DELTA
    mod4=(T0*T0+DELTA*DELTA)**2           # |a|^4
    beta=-2*a2r/mod4                       # -2 Re(1/a^2)
    kappa=1/mod4                           # |a|^-4
    a=[]
    for m in range(Mmax):
        v=gam[m]
        if m>=1: v+=-beta*(2*m)*(2*m-1)*gam[m-1]
        if m>=2: v+= kappa*(2*m)*(2*m-1)*(2*m-2)*(2*m-3)*gam[m-2]
        a.append(v)
    return a
def jensen_maximag(a,d,n):
    co=[mp.binomial(d,k)*a[n+k] for k in range(d,-1,-1)]
    sc=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/sc for c in co]
    rts=mp.polyroots(co,maxsteps=400,extraprec=250)
    return max(abs(mp.im(r)) for r in rts)/max(abs(r) for r in rts)
def first_failure(a,dmax):
    # returns (d*, n*, maximag) or None
    for d in range(2,dmax+1):
        worst=mp.mpf(0); wn=None
        for n in range(0,Mmax-d):
            v=jensen_maximag(a,d,n)
            if v>worst: worst,wn=v,n
            if v>THR:
                return d,n,v,worst  # first failing d, at shift n
    return None
print(f"d*(T0) experiment  DELTA={DELTA}  moments to m={Mmax-1}  thr={mp.nstr(THR,2)}",flush=True)
print(" T0 :  d*   n*    max|Im|@d*     (d*-1 fully hyperbolic below thr)",flush=True)
res=[]
for T0 in [3,4,5,6,8,10]:
    a=target_moments(T0)
    dmax=min(Mmax-1, int(3*T0*T0)+5)
    f=first_failure(a,dmax)
    if f is None:
        print(f" {T0:3d}:  none up to d={dmax} (moment-limited)",flush=True)
        continue
    dstar,nstar,vfail,_=f
    # two-sided: confirm d*-1 all hyperbolic
    below_ok=True; below_worst=mp.mpf(0)
    if dstar>2:
        for n in range(0,Mmax-(dstar-1)):
            vv=jensen_maximag(a,dstar-1,n); below_worst=max(below_worst,vv)
            if vv>THR: below_ok=False;break
    res.append((mp.mpf(T0),dstar))
    print(f" {T0:3d}:  {dstar:3d}  {nstar:3d}   {mp.nstr(vfail,4)}   d*-1 max|Im|={mp.nstr(below_worst,3)} {'OK(<thr)' if below_ok else 'FAIL-not-clean'}",flush=True)
print("",flush=True)
if len(res)>=3:
    print("fit log(d*) = A + B log(T0):",flush=True)
    xs=[mp.log(t) for t,_ in res]; ys=[mp.log(mp.mpf(d)) for _,d in res]
    n=len(xs); sx=sum(xs); sy=sum(ys); sxx=sum(x*x for x in xs); sxy=sum(x*y for x,y in zip(xs,ys))
    B=(n*sxy-sx*sy)/(n*sxx-sx*sx); A=(sy-B*sx)/n
    print(f"  exponent B = {mp.nstr(B,5)}  (NULL: B=2.0)  intercept A={mp.nstr(A,4)} c=e^A={mp.nstr(mp.e**A,4)}",flush=True)
    print("  residuals log d* - fit:",flush=True)
    for (t,d),x,y in zip(res,xs,ys):
        print(f"    T0={mp.nstr(t,3)}: d*={d}  pred={mp.nstr(mp.e**(A+B*x),4)}  resid={mp.nstr(y-(A+B*x),3)}",flush=True)
print("DONE",flush=True)

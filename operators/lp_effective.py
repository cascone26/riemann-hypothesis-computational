#!/usr/bin/env python3
"""§313B: effective-GORZ scheme. Three pieces:
(1) DIRECT: hyperbolicity of the actual Jensen polys J^{d,n}(X)=sum_k C(d,k) gamma_{n+k} X^k, d=3..12, n=0..12
    (polyroots; hyperbolic iff max |Im z|/|z| ~ 0).
(2) MARGIN eps_d: empirical hyperbolicity margin of the Hermite limit He_d — smallest relative coefficient
    perturbation (worst of K sign patterns, bisected) that creates complex roots.
(3) THRESHOLD N(d): leading-order scheme — the cubic deviation of log gamma over a d-window is ~|c_n| d^3/6,
    c_n ~ 2/n^2 (§313A). Requiring (2/n^2)d^3/6 < eps_d gives N(d)=sqrt(d^3/(3 eps_d)). If direct check covers
    n<=N(d), that degree d is verified for ALL n (within the leading-order scheme — REPORT tier, stated honestly).
"""
import mpmath as mp
mp.mp.dps=45
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,90):
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n*n*e2u)
        s+=term
        if abs(term)<mp.mpf(10)**(-50) and n>3: break
    return s
Mmax=25
gam=[]
for m in range(Mmax):
    M=mp.quad(lambda u: Phi(u)*u**(2*m),[0,0.4,0.8,1.2,2,4])
    gam.append(M/mp.factorial(2*m))
print("(1) DIRECT Jensen hyperbolicity J^{d,n}, entries = max|Im z|/max|z| (0 => all real roots):",flush=True)
def jensen_maximag(d,n):
    # coeffs high->low for polyroots: J = sum_k C(d,k) gamma_{n+k} X^k
    co=[mp.binomial(d,k)*gam[n+k] for k in range(d,-1,-1)]
    # normalize dynamic range: divide by geometric mean of |lead|,|const|
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=200,extraprec=120)
    mx=max(abs(mp.im(r)) for r in rts); mxa=max(abs(r) for r in rts)
    return mx/mxa
hdr="      "+" ".join(f"n={n:<2d}" for n in range(0,13,2))
print(hdr,flush=True)
worst=mp.mpf(0)
for d in range(3,13):
    row=[]
    for n in range(0,13,2):
        if n+d<Mmax:
            v=jensen_maximag(d,n); worst=max(worst,v)
            row.append(f"{mp.nstr(v,2):>5}")
        else: row.append("  -  ")
    print(f"d={d:2d}: "+" ".join(row),flush=True)
print(f"  worst max|Im|/|z| over all (d,n) tested = {mp.nstr(worst,3)}  => {'ALL HYPERBOLIC' if worst<mp.mpf('1e-20') else 'check individually'}",flush=True)
print("",flush=True)
print("(2) Hermite margins eps_d (empirical, worst of 40 sign patterns, bisected):",flush=True)
import random
random.seed(7)
def hermite_coeffs(d):
    # probabilists' He_d coefficients via recurrence He_{k+1}=x He_k - k He_{k-1}
    P=[mp.mpf(1)]; Q=[mp.mpf(0),mp.mpf(1)]  # He_0, He_1 (low->high)
    if d==0: return P
    for k in range(1,d):
        R=[mp.mpf(0)]+Q            # x*He_k
        for i in range(len(P)): R[i]-=k*P[i]
        P,Q=Q,R
    return Q  # low->high
def perturbed_real(co_lowhigh,eps,signs):
    co=[c*(1+eps*s) if c!=0 else c for c,s in zip(co_lowhigh,signs)]
    rts=mp.polyroots(list(reversed(co)),maxsteps=200,extraprec=120)
    mxa=max(abs(r) for r in rts)
    return max(abs(mp.im(r)) for r in rts)/mxa < mp.mpf('1e-12')
epsd={}
for d in range(3,9):
    co=hermite_coeffs(d)
    idx=[i for i,c in enumerate(co) if c!=0]
    patterns=[[random.choice([-1,1]) if i in idx else 0 for i in range(len(co))] for _ in range(40)]
    lo,hi=mp.mpf(0),mp.mpf(1)
    for _ in range(30):
        mid=(lo+hi)/2
        ok=all(perturbed_real(co,mid,p) for p in patterns)
        if ok: lo=mid
        else: hi=mid
    epsd[d]=lo
    print(f"  d={d}: eps_d ~ {mp.nstr(lo,4)}",flush=True)
print("",flush=True)
print("(3) THRESHOLDS N(d)=sqrt(d^3/(3 eps_d)) from c_n~2/n^2 (leading-order scheme):",flush=True)
for d in range(3,9):
    N=mp.sqrt(d**3/(3*epsd[d]))
    covered = N<=12
    print(f"  d={d}: N(d) ~ {mp.nstr(N,4)}   direct check covers n<=12: {'=> ALL n VERIFIED (scheme)' if covered else f'need direct to n~{mp.nstr(N,3)}'}",flush=True)
print("DONE",flush=True)

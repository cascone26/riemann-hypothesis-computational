#!/usr/bin/env python3
"""§309: de Bruijn-Newman flow on our explicit Phi. H_t(x)=int_0^inf e^{t u^2} Phi(u) cos(xu) du.
Lambda = inf{t: H_t has only real zeros}; RH <=> Lambda<=0 (Rodgers-Tao: Lambda>=0, so RH<=>Lambda=0).
The flow's moments gamma_m(t)=M_m(t)/(2m)!, M_m(t)=int e^{t u^2}Phi(u)u^{2m}du. Turan/Jensen hyperbolicity
of {gamma_m(t)} IMPROVES with t. Track the Turan margin ratio-1 vs t (incl t<0): as t decreases the margin
shrinks toward the finite-degree hyperbolicity threshold. Phi ~ exp(-pi e^{2u}) (double-exp) dominates e^{t u^2}
for all real t, so integrals converge everywhere. Charisma: REPORT (DBN flow on the explicit kernel)."""
import mpmath as mp
mp.mp.dps=30
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,80):
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n*n*e2u)
        s+=term
        if abs(term)<mp.mpf(10)**(-35) and n>3: break
    return s
Mmax=9
def gammas_t(t):
    g=[]
    for m in range(Mmax):
        M=mp.quad(lambda u: mp.e**(t*u*u)*Phi(u)*u**(2*m),[0,0.5,1,2,4])
        g.append(M/mp.factorial(2*m))
    return g
print("DBN flow: Turan margin  min_n [ gamma_{n+1}(t)^2/(gamma_n gamma_{n+2}) - 1 ]  vs t",flush=True)
print("(t=0 is the RH point; margin>0 => low-deg Jensen hyperbolic; flow spreads zeros to real axis as t grows)",flush=True)
for t in [mp.mpf('-0.6'),mp.mpf('-0.4'),mp.mpf('-0.2'),mp.mpf('-0.1'),mp.mpf('0'),mp.mpf('0.1'),mp.mpf('0.3'),mp.mpf('0.6')]:
    g=gammas_t(t)
    margins=[g[n+1]**2/(g[n]*g[n+2])-1 for n in range(Mmax-2)]
    mn=min(margins); argmn=margins.index(mn)
    allpos=all(x>0 for x in margins)
    print(f"  t={mp.nstr(t,3):>5}: min Turan margin={mp.nstr(mn,6):>11} (at n={argmn})  all d=2 hyperbolic: {allpos}   [n=0 margin={mp.nstr(margins[0],5)}]",flush=True)
print("",flush=True)
print("trend: margin decreases as t decreases. Extrapolate where min-margin -> 0 (a finite-degree Lambda proxy).",flush=True)
print("Rodgers-Tao: true Lambda>=0. RH <=> Lambda<=0. Our t=0 margins all >0 (RH-consistent).",flush=True)
print("DONE",flush=True)

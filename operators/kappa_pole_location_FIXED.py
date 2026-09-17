#!/usr/bin/env python3
"""§973 (CORRECTED) κ = #pole-locations at s=1, cross-Selberg — proper primitive Kronecker characters.

Supersedes kappa_equals_pole_order.py, whose legendre_char used jacobi(n,q) — INVALID for even modulus
(q=4,8): the Jacobi symbol requires an odd denominator, so χ_4/χ_8 there were garbage and spuriously gave
κ=1 (the workflow verifier correctly caught this). Here the real primitive quadratic characters are the
Kronecker symbols: χ_4=(-1/n), χ_8=(2/n) even, χ_8'=(-2/n) odd, and (·/q) for odd q. Result (VERIFIED):
ζ→κ=1 (pole), every entire χ→κ=0. HANDLE."""
import numpy as np
from scipy.special import digamma
def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=8000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMv=LAM[2:]
def chi4(n):    n%=4; return {1:1,3:-1}.get(n,0)              # (-1/n) odd
def chi8_2(n):  n%=8; return {1:1,7:1,3:-1,5:-1}.get(n,0)     # (2/n) even
def chi8_m2(n): n%=8; return {1:1,3:1,5:-1,7:-1}.get(n,0)     # (-2/n) odd
def jacobi_odd(a,n):
    a%=n; r=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): r=-r
        a,n=n,a
        if a%4==3 and n%4==3: r=-r
        a%=n
    return r if n==1 else 0
def chi_odd(q): return lambda n: jacobi_odd(n%q,q)            # (·/q) for ODD q only
def weil_kappa(chi,q,a_par,fourW,cen_lo=14,cen_hi=44,frac=0.8,dr=0.004,is_zeta=False):
    W=fourW/4; L=cen_hi-cen_lo; N=max(4,round(frac*L*2*W/np.pi))
    RG=np.arange(-(cen_hi+120),(cen_hi+120),dr)
    OM=np.real(digamma((a_par+0.5)/2+0.5j*RG))+np.log(q/np.pi)
    CEN=np.linspace(cen_lo,cen_hi,N)
    def fej(x):
        y=W*x;o=np.ones_like(y);nz=np.abs(y)>1e-12;o[nz]=(np.sin(y[nz])/y[nz])**2;return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); G=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi); mask=LOGN<fourW; ln=LOGN[mask]
    chiv=np.array([1.0 if is_zeta else chi(int(n)) for n in NS[mask]])
    wts=LAMv[mask]*chiv/np.sqrt(NS[mask])
    cosu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*wts[None,:]).sum(axis=1)
    M=arch-0.5*(prime+prime.T)
    gw,gv=np.linalg.eigh(G); k=gw>gw.max()*1e-12; B=gv[:,k]/np.sqrt(gw[k])
    ev=np.linalg.eigvalsh(B.T@M@B); s=max(1,abs(ev).max()); return int((ev<-1e-9*s).sum())
if __name__=="__main__":
    print(f"{'L':>12} {'q':>3} {'a':>2} {'pole':>4} | κ@4W=5 κ@4W=6")
    for name,chi,q,a,isz in [("zeta",None,1,0,True),("chi4(odd)",chi4,4,1,False),
            ("chi8_2(even)",chi8_2,8,0,False),("chi8_m2(odd)",chi8_m2,8,1,False),
            ("chi3(odd)",chi_odd(3),3,1,False),("chi5(even)",chi_odd(5),5,0,False),
            ("chi7(odd)",chi_odd(7),7,1,False),("chi11(odd)",chi_odd(11),11,1,False)]:
        ks=[weil_kappa(chi,q,a,fw,is_zeta=isz) for fw in (5.0,6.0)]
        pole=1 if isz else 0
        print(f"{name:>12} {q:>3} {a:>2} {pole:>4} |  {ks[0]:>4}   {ks[1]:>4}  {'OK' if all(k==pole for k in ks) else '*** MISMATCH'}")

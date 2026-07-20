#!/usr/bin/env python3
"""
critical_eigvec_id.py -- ATTACK A2: identify the CRITICAL (near-kernel) eigenvector of tau
in closed form.

The min-|eigenvalue| eigenvector g_min is the TIGHTEST direction of the positivity wall
(where A+B-P is closest to failing). weil_cancellation_analysis showed it converges to a
fixed limit function as N->inf. If we can NAME g_lim (constant? Lorentzian 1/(n^2+a)?
exponential? prolate spheroidal?), we can try to prove B(g_lim) >= |A(g_lim)|+P(g_lim)
analytically -> closes the wall on the critical direction. This is a proof attack, not a stat.

Build tau (Groskin cutoff-free closed form, float64 -- eigenvector SHAPE needs speed not 60
digits), get the near-kernel eigenvector, and:
  1. symmetrize to even part (the construction's parity sector); normalize.
  2. examine coefficient decay vs n: fit to candidate laws (const, e^{-a|n|}, 1/(n^2+a^2),
     sech, gaussian) -> which matches?
  3. check CONVERGENCE of the coefficients across N (a real limit function must stabilize).
"""
import numpy as np
from scipy.special import digamma
import mpmath as _mp
def polygamma1(z):  # complex trigamma via mpmath (scipy's can't take complex)
    v=_mp.polygamma(1, complex(z)); return complex(v)
import json, sys

TWO_PI = 2*np.pi

def prime_powers_up_to(c):
    primes,x=[],2
    while x<=c:
        if all(x%p for p in primes): primes.append(x)
        x+=1
    out=[]
    for p in primes:
        q=p
        while q<=c: out.append((q,p)); q*=p
    return out

def geom_sums(n,L):
    w=2*np.pi*n/L; w2=w*w
    ks=np.arange(0,80); ck=2*ks+0.5; e=np.exp(-ck*L); den=ck*ck+w2
    gS=np.sum(e/den); gCC=np.sum(e*w2/(ck*den)) if n!=0 else 0.0
    gX1=np.sum(e*ck/den); gX2=np.sum(e*(ck*ck-w2)/(den*den))
    return gS,gCC,gX1,gX2

def build_tau_np(N,c):
    L=np.log(c); PI=np.pi; psi_q=digamma(0.25)
    S=np.zeros(N+1);CC=np.zeros(N+1);XC=np.zeros(N+1)
    for n in range(N+1):
        w=2*PI*n/L; z=0.25+1j*PI*n/L
        psi=digamma(z); psi1=polygamma1(z)
        gS,gCC,gX1,gX2=geom_sums(n,L)
        if n==0: S[n]=0;CC[n]=0
        else: S[n]=psi.imag/2-w*gS; CC[n]=-(psi.real-psi_q)/2+gCC
        XC[n]=psi1.real/4-L*gX1-gX2
    sp2=16*PI*PI;l2=L*L; pref02=32*L*np.sinh(L/4)**2
    eL=np.exp(L);kappa=np.log(4*PI*(eL-1)/(eL+1))+np.euler_gamma
    U=np.exp(L/2);J=-2*np.log(U+1)+np.log(U*U+1)+2*np.arctan(U)+np.log(2)-PI/2
    pdata=prime_powers_up_to(c)
    wts=np.array([np.log(p)*q**-0.5 for (q,p) in pdata]); pos=np.array([np.log(q) for (q,p) in pdata])
    def Ssg(nn): return S[nn] if nn>=0 else -S[-nn]
    D=2*N+1; tau=np.zeros((D,D))
    ns=np.arange(-N,N+1)
    for i in range(D):
        n=ns[i]
        for j in range(i,D):
            m=ns[j]
            num=l2-sp2*m*n;den=(l2+sp2*m*m)*(l2+sp2*n*n);W02=pref02*num/den
            if n==m: WR=kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]
            else: WR=(Ssg(m)-Ssg(n))/(PI*(n-m))
            if n==m: qv=2*(1-pos/L)*np.cos(2*PI*n*pos/L)
            else: qv=(np.sin(2*PI*m*pos/L)-np.sin(2*PI*n*pos/L))/(PI*(n-m))
            Wp=np.sum(wts*qv)
            v=W02-WR-Wp; tau[i,j]=v;tau[j,i]=v
    return tau,ns

def candidate_fits(nn, a):
    """a = |coeff| at |n|=nn (nn>=0), normalized a[0]=1. Fit decay laws, return best."""
    x=np.arange(len(a)); y=a/a[0]
    mask=y>1e-6
    x,y=x[mask],y[mask]
    fits={}
    # exponential e^{-b n}
    if len(x)>3:
        b=np.polyfit(x,np.log(y),1); fits["exp e^{-b n}"]=(-b[0], np.exp(np.polyval(b,x)))
    # lorentzian 1/(n^2+a^2): log y = -log(n^2+a^2)+const -> fit a^2
    try:
        from scipy.optimize import curve_fit
        f_lor=lambda n,A,a2: A/(n*n+a2)
        popt,_=curve_fit(f_lor,x,y,p0=[y[0],1.0],maxfev=5000); fits["lorentz 1/(n^2+a^2)"]=(popt[1],f_lor(x,*popt))
        f_gau=lambda n,A,s: A*np.exp(-(n/s)**2)
        popt2,_=curve_fit(f_gau,x,y,p0=[y[0],2.0],maxfev=5000); fits["gauss e^{-(n/s)^2}"]=(popt2[1],f_gau(x,*popt2))
        f_sech=lambda n,A,s: A/np.cosh(n/s)
        popt3,_=curve_fit(f_sech,x,y,p0=[y[0],2.0],maxfev=5000); fits["sech(n/s)"]=(popt3[1],f_sech(x,*popt3))
    except Exception as e:
        fits["_err"]=str(e)
    # report RMS of each on log scale
    scores={}
    for k,(param,pred) in fits.items():
        if k=="_err": continue
        r=np.sqrt(np.mean((np.log(np.maximum(pred,1e-12))-np.log(y))**2))
        scores[k]={"param":float(param),"logRMS":float(r)}
    return scores

if __name__=="__main__":
    c=int(sys.argv[1]) if len(sys.argv)>1 else 14
    Ns=[int(x) for x in sys.argv[2].split(",")] if len(sys.argv)>2 else [48,96,160,256]
    out={"c":c,"runs":{}}
    coeff_by_N={}
    for N in Ns:
        tau,ns=build_tau_np(N,c)
        ev,evec=np.linalg.eigh(tau)
        k0=int(np.argmin(np.abs(ev)))
        v=evec[:,k0]; D=2*N+1
        even=np.array([(v[i]+v[D-1-i])/2 for i in range(D)])
        # even part indexed by |n|=0..N
        a=np.array([abs(even[N+k]) for k in range(N+1)])
        if a[0]==0: a[0]=max(a.max(),1e-12)
        a=a/a[0]
        coeff_by_N[N]=a
        scores=candidate_fits(np.arange(N+1), a)
        best=min(scores.items(), key=lambda kv: kv[1]["logRMS"]) if scores else ("none",{})
        print(f"\n=== N={N}: near-kernel eigval={ev[k0]:.2e}, |eigval|rank0 ===", flush=True)
        print(f"  coeff |c_n| (n=0..8): {[f'{x:.3f}' for x in a[:9]]}")
        for k,s in sorted(scores.items(), key=lambda kv: kv[1]['logRMS']):
            print(f"    {k:24s} param={s['param']:.4f}  logRMS={s['logRMS']:.3f}")
        print(f"  BEST FIT: {best[0]}")
        out["runs"][N]={"near_eigval":float(ev[k0]),"coeff_0_8":a[:9].tolist(),
                        "fits":scores,"best":best[0]}
    # convergence of coefficients across N (compare first 9 coeffs)
    print("\n"+"="*70); print("COEFFICIENT CONVERGENCE across N (n=0..8):")
    for k in range(9):
        row=[f"{coeff_by_N[N][k]:.4f}" for N in Ns]
        print(f"  c_{k}: "+"  ".join(f"N{N}={coeff_by_N[N][k]:.4f}" for N in Ns))
    p="/Users/scones/projects/riemann/results/critical_eigvec_id.json"
    json.dump(out,open(p,"w"),indent=2)
    print(f"\nsaved -> {p}")

#!/usr/bin/env python3
"""§310: structural handle on the LP/Jensen front (NOT numerics-for-their-own-sake): the RATE the Jensen polys
approach Hermite (= high-degree hyperbolicity margin) is set by Phi's modular e^{2u} SADDLE. Compute gamma_m to
m=20, extract the Turan log-ratio r_m=log(gamma_{m-1}gamma_{m+1}/gamma_m^2) and its decay, and DERIVE the rate from
the saddle point u*(m) of M_m=int Phi u^{2m}. If exact r_m matches the saddle prediction, we have an ANALYTIC handle
on the margin (from the modular exponent), not a fit. Charisma: REPORT + a saddle-vs-exact check."""
import mpmath as mp
mp.mp.dps=50
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,90):
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n*n*e2u)
        s+=term
        if abs(term)<mp.mpf(10)**(-55) and n>3: break
    return s
Mmax=21
gam=[]
for m in range(Mmax):
    M=mp.quad(lambda u: Phi(u)*u**(2*m),[0,0.4,0.8,1.2,2,4])
    gam.append(M/mp.factorial(2*m))
print("Turan log-ratio r_m = log(gamma_{m-1} gamma_{m+1}/gamma_m^2)  (>0 => log-concave/hyperbolic d=2), and decay:",flush=True)
prev=None
for m in range(1,Mmax-1):
    r=mp.log(gam[m-1]*gam[m+1]/gam[m]**2)   # note: this is NEGATIVE of Turan margin sign; r<0 => log-concave
    # actually log-concave <=> gamma_m^2>=gamma_{m-1}gamma_{m+1} <=> ratio(m-1,m+1)/m^2 <=1 <=> r<=0
    rate = None if prev is None else prev/r
    print(f"  m={m:2d}: r_m={mp.nstr(r,7):>12}  (<0=hyperbolic)  r_{{m-1}}/r_m={mp.nstr(rate,5) if rate else '-'}  m^2*|r_m|={mp.nstr(m*m*abs(r),6)}",flush=True)
    prev=r
print("",flush=True)
# Saddle-point prediction: M_m=int Phi(u)u^{2m}du, dominant Phi~2pi^2 e^{9u/2} e^{-pi e^{2u}} (n=1, large u).
# f(u)=2m ln u + 9u/2 - pi e^{2u}; saddle 2m/u + 9/2 = 2 pi e^{2u}. Solve u*(m), Laplace approx gamma_m, compare.
print("saddle u*(m) and Laplace gamma_m vs exact (validates modular-saddle control of the margin):",flush=True)
def usaddle(m):
    return mp.findroot(lambda u: 2*m/u+mp.mpf(9)/2-2*mp.pi*mp.e**(2*u), mp.mpf(0.5)+mp.log(m+1)/2)
for m in [6,10,14,18]:
    us=usaddle(m)
    fpp=-2*m/us**2-4*mp.pi*mp.e**(2*us)
    Mlap=2*mp.pi**2*mp.e**(mp.mpf(9)*us/2)*mp.e**(-mp.pi*mp.e**(2*us))*us**(2*m)*mp.sqrt(2*mp.pi/(-fpp))
    glap=Mlap/mp.factorial(2*m)
    print(f"  m={m}: u*={mp.nstr(us,6)}  gamma_saddle={mp.nstr(glap,6)}  gamma_exact={mp.nstr(gam[m],6)}  ratio={mp.nstr(glap/gam[m],6)}",flush=True)
print("DONE",flush=True)

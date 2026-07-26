#!/usr/bin/env python3
"""§307: FRESH family (leaves the positivity-criterion wall) — Jensen polynomials / GORZ 2019.
Xi(t)=int_0^inf Phi(u)cos(tu)du = sum_m (-1)^m gamma_m t^{2m}, gamma_m = M_m/(2m)! > 0 (Phi>0),
M_m = int_0^inf Phi(u) u^{2m} du.  Jensen poly J^{d,n}(X)=sum_{k=0}^d C(d,k) gamma_{n+k} X^k.
RH <=> every J^{d,n} is HYPERBOLIC (all real roots) [Polya; GORZ proved it for fixed d, all large n].
d=2 hyperbolic <=> gamma_{n+1}^2 >= gamma_n gamma_{n+2} (Turan / log-concavity). Compute gamma_m, check
Turan (d=2) and d=3 hyperbolicity. Charisma: REPORT (numerics on a mainline RH-equivalent)."""
import mpmath as mp
mp.mp.dps=40
def Phi(u):
    s=mp.mpf(0)
    e2u=mp.e**(2*u)
    for n in range(1,60):
        t=mp.pi*n*n*e2u
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-t)
        s+=term
        if abs(term)<mp.mpf(10)**(-45) and n>3: break
    return s
Mmax=13
print("computing moments M_m = int_0^inf Phi(u)u^{2m}du (Phi>0, Polya kernel):",flush=True)
gam=[]
for m in range(Mmax):
    M=mp.quad(lambda u: Phi(u)*u**(2*m),[0,0.5,1,2,4])
    g=M/mp.factorial(2*m)
    gam.append(g)
    print(f"  gamma_{m} = {mp.nstr(g,12)}  (>0 ? {g>0})",flush=True)
print("",flush=True)
# d=2 Jensen = Turan / log-concavity: gamma_{n+1}^2 >= gamma_n gamma_{n+2}
print("d=2 Jensen (Turan): gamma_{n+1}^2 - gamma_n*gamma_{n+2}  (>=0 => hyperbolic, RH-consistent):",flush=True)
allpos=True
for n in range(Mmax-2):
    t=gam[n+1]**2-gam[n]*gam[n+2]
    ratio=gam[n+1]**2/(gam[n]*gam[n+2])
    ok=t>=0
    allpos=allpos and ok
    print(f"  n={n}: Turan={mp.nstr(t,6)}  ratio(>=1)={mp.nstr(ratio,8)}  {'OK' if ok else 'VIOLATION!'}",flush=True)
print(f"  => all d=2 Jensen hyperbolic in range: {allpos}",flush=True)
# d=3 Jensen: J = gamma_n + 3 gamma_{n+1}X + 3 gamma_{n+2}X^2 + gamma_{n+3}X^3, check 3 real roots via discriminant>0
print("d=3 Jensen hyperbolicity (cubic discriminant > 0 => 3 real roots):",flush=True)
for n in range(Mmax-3):
    a=gam[n+3]; b=3*gam[n+2]; c=3*gam[n+1]; d=gam[n]   # a X^3 + b X^2 + c X + d  (coeffs high->low)
    disc=18*a*b*c*d-4*b**3*d+b*b*c*c-4*a*c**3-27*a*a*d*d
    print(f"  n={n}: disc={mp.nstr(disc,6)}  {'HYPERBOLIC (3 real roots)' if disc>0 else 'NOT (complex roots) -- RH VIOLATION' }",flush=True)
print("DONE",flush=True)

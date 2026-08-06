"""Jensen hyperbolicity with the CORRECT sequence: G(x)=Sum_m [M(m)/(2m)!] x^m, RH <=> G in
Laguerre-Polya <=> Jensen polys of alpha(m)=M(m)*m!/(2m)! all hyperbolic. High-precision roots."""
import mpmath as mp
from math import comb
mp.mp.dps = 100
def Phi(u):
    u=mp.mpf(u); tot=mp.mpf(0)
    e4u=mp.e**(4*u); e9u=mp.e**(9*u); e5u=mp.e**(5*u)
    for n in range(1,100):
        term=(2*mp.pi**2*n**4*e9u-3*mp.pi*n**2*e5u)*mp.e**(-mp.pi*n**2*e4u)
        tot+=term
        if abs(term)<mp.mpf(10)**(-110) and n>5: break
    return tot
NM=24
M=[mp.quad(lambda u: u**(2*n)*Phi(u), [0,0.5,1,2,4]) for n in range(NM)]
alpha=[M[m]*mp.factorial(m)/mp.factorial(2*m) for m in range(NM)]
print("alpha(0..5):",[mp.nstr(alpha[m],5) for m in range(6)], " all>0:",all(a>0 for a in alpha))
def check(d,n):
    coeffs=[mp.mpf(comb(d,j))*alpha[n+j] for j in range(d+1)]
    roots=mp.polyroots(coeffs[::-1], maxsteps=300, extraprec=600)
    mi=max(abs(mp.im(r)) for r in roots); mr=max(abs(mp.re(r)) for r in roots)
    negreal=all(mp.re(r)<0 for r in roots)  # LP: all negative real
    return float(mi/mr), negreal
print(f"\n{'d':>3}{'n':>3}  {'max|Im|/|Re|':>14}  {'all neg-real?':>13}  hyperbolic?")
allhyp=True
for d in [3,5,8,10,12,15,18,20]:
    for n in [0,3]:
        if n+d>=NM: continue
        ratio,neg=check(d,n); hyp=ratio<1e-30
        allhyp=allhyp and hyp
        print(f"{d:>3}{n:>3}  {ratio:>14.2e}  {str(neg):>13}  {'YES' if hyp else 'NO'}")
print(f"\nALL hyperbolic (RH-consistent, GORZ reproduced): {allhyp}")

# GORZ Hermite limit: appropriately centered+scaled J^{d,n} roots -> Hermite H_d roots as n->inf.
# Use the natural centering (from log-derivative of alpha) and unit-variance scaling.
import numpy as np
from numpy.polynomial.hermite import hermroots
print("\n=== GORZ Hermite limit check (d=5): centered+scaled Jensen roots -> Hermite H_5 roots ===")
hd=np.array(sorted(hermroots([0]*5+[1])))     # physicists' H_5 roots
hd_n=hd/np.std(hd)                             # unit-variance
d=5
for n in [2,6,10,14,18]:
    if n+d>=NM: continue
    coeffs=[mp.mpf(comb(d,j))*alpha[n+j] for j in range(d+1)]
    roots=sorted(float(mp.re(r)) for r in mp.polyroots(coeffs[::-1],maxsteps=300,extraprec=600))
    r=np.array(roots); rc=(r-r.mean()); rc=rc/rc.std()   # center + unit variance
    err=np.max(np.abs(rc-hd_n))
    print(f"  n={n:2}: max|scaled Jensen root - scaled Hermite root| = {err:.4f}")
print("  (-> 0 as n grows confirms the GORZ Hermite universality; validates setup)")
